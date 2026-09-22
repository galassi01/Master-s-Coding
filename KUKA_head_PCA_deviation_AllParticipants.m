function KUKA_head_PCA_deviation_AllParticipants(participant_nums, base_dir, filt, opts)
%% ========================================================================
%  HEAD-VS-FLANGE PCA ROTATION-AXIS DEVIATION -- ALL-PARTICIPANT BATCH
% ========================================================================
%  Complements (does NOT replace) KUKA_head_kinematics_summary_
%  AllParticipants.m's magnitude-based off-axis metrics (OffAxis_Resultant_
%  RMS/CV, OffAxis_to_Task_RMS_Ratio). Those are all MAGNITUDE measures --
%  they scale with how much movement happened, which is exactly why
%  frequency/amplitude confound them. This script instead asks a
%  SCALE-INVARIANT geometric question: irrespective of how much the head
%  moved, WHICH AXIS did it move around, and how far off is that axis from
%  the axis the flange actually delivered that trial.
%
%  METHOD
%  ------
%  For each (Participant, Direction, ParameterSet, Condition) group:
%    1. Identify every trimmed_*.mat file belonging to that group. Each
%       PHYSICAL trial is written out once per Muscle bank (1 and 2) with
%       IDENTICAL head/flange IMU data (only the EMG channels differ), so
%       duplicates are collapsed to one file per physical trial (lowest
%       Muscle number kept) before pooling -- otherwise every physical
%       trial's cycles would be double-counted.
%    2. POOL (concatenate) the CYCLE-TRIMMED, resampled head and flange
%       gyro samples (data.imu_head.six_axis_resamp(:,4:6) / data.imu_
%       flange.six_axis_resamp(:,4:6) -- the same "_Resamp" block the
%       off-axis magnitude script uses) across every physical trial in
%       the group. A single trial may only contain a handful of cycles;
%       pooling repeat trials (e.g. 2x20 or 4x10 cycles) gets a stable
%       ~10-20+ cycle sample for the PCA without averaging anything, and
%       naturally weights each trial by how many samples/cycles it
%       actually contributed -- unlike averaging each trial's own PC1
%       vector, which is NOT done here (eigenvectors are sign-ambiguous,
%       so naively averaging them can partially cancel, and it implicitly
%       weights every trial equally regardless of cycle count -- both
%       wrong). One PCA is run on the pooled sample instead.
%    3. Run 3D PCA (eig of the covariance, mean-centered) separately on
%       the pooled head gyro and the pooled flange gyro -> head_PC1 and
%       flange_PC1, each a unit vector in (Gx,Gy,Gz) space, plus VAF
%       (fraction of variance explained by PC1 -- a reliability/
%       confidence indicator: VAF near 1 means a well-defined dominant
%       rotation axis; VAF near 1/3 means the motion is nearly isotropic
%       and PC1's direction is not meaningfully defined).
%    4. DeviationAngle_3D_deg = acosd(|dot(head_PC1, flange_PC1)|) -- the
%       angle (0-90 deg) between the head's actual dominant rotation axis
%       and the axis the flange ACTUALLY delivered that trial (not the
%       idealized/nominal design axis -- this is what was really
%       delivered, imperfections included). The absolute value removes
%       the eigenvector sign ambiguity (PC1 is only defined up to +/-1),
%       so no reference-vector sign-fix is needed for this number.
%    5. IN-PLANE vs OUT-OF-PLANE SPLIT. The flange is mechanically
%       confined to the horizontal roll/pitch (Gx-Gy) plane by design --
%       Gz is never a task axis for ML/AP/diagonal trials (same
%       assumption already used throughout the kinematics scripts; this
%       run's Flange_OutOfPlane_Angle_deg is a free QA check on that
%       assumption -- it should come out near 0). Splitting head_PC1
%       against that plane separates two mechanistically different
%       errors that a single 3D angle would blur together:
%         OutOfPlane_Angle_deg        = asind(|head_PC1_z|) -- how much
%           of the head's dominant axis left the roll/pitch plane
%           entirely (yaw recruitment), regardless of bearing.
%         InPlane_Azimuth_Angle_deg   = angle, WITHIN that plane, between
%           the head's projected axis and the direction's nominal in-
%           plane axis (unsigned, 0-90 deg) -- "given it stayed in the
%           right plane, how far off was the bearing."
%         InPlane_Azimuth_SignedDeg   = same, signed (-90 to +90 deg) --
%           requires sign-fixing head_PC1 first (dot against the nominal
%           axis, flip if negative) since a signed angle on a sign-
%           ambiguous vector is meaningless; useful for seeing whether
%           deviation is consistently biased to one side rather than
%           random. NaN when the head's motion has (numerically) no
%           in-plane component at all.
%       Nominal per-direction in-plane axis (same convention already
%       used for diagonal-trial sign-fixing elsewhere in this pipeline):
%         ML   -> [1, 0, 0]              AP   -> [0, 1, 0]
%         DiaR -> [1, 1, 0]/sqrt(2)      DiaL -> [1,-1, 0]/sqrt(2)
%
%    6. PER-CALIBRATED-AXIS DEVIATION VECTOR (DeviationVec_Gx/Gy/Gz).
%       head_PC1/flange_PC1 are already expressed in the calibrated
%       (Gx,Gy,Gz) frame from the two-pose calibration (createCalibration
%       in KUKA_analysis_prep_consume_ref.m) -- PCA doesn't change
%       coordinate systems, it just picks the empirically dominant
%       direction WITHIN that frame, which is deliberately preferred over
%       assuming the diagonals are exactly a 45/-45 deg split of Gx/Gy
%       (real mechanical mounting/calibration is never exactly 45 deg).
%       To get a per-axis deviation, head_PC1 is sign-fixed against
%       flange_PC1 itself (NOT the nominal design axis used for
%       OutOfPlane/InPlane above -- those need the nominal reference
%       because they're defined relative to the idealized direction, but
%       a head-vs-flange deviation vector needs both sides referenced to
%       what the flange ACTUALLY delivered, or the component-wise
%       difference would be meaningless):
%         DeviationVec_Gx = head_PC1_flangealigned_x - flange_PC1_x
%         DeviationVec_Gy = head_PC1_flangealigned_y - flange_PC1_y
%         DeviationVec_Gz = head_PC1_flangealigned_z - flange_PC1_z
%       These are unitless (both inputs are unit vectors) and tell you
%       WHICH calibrated axis the head's rotation leaked into relative to
%       the flange, and which direction (sign) -- e.g. a diagonal trial
%       might show +Gx, -Gy, +Gz rather than assuming any fixed split.
%       DeviationAngle_3D_deg (above) is the norm/magnitude of essentially
%       this same comparison collapsed to one number; this vector is the
%       axis-resolved version of the same underlying comparison.
%
%  OUTPUT LEVELS (mirrors KUKA_head_kinematics_summary_AllParticipants.m)
%    Level A  Head_PCA_Deviation_PerParticipant.csv
%             One row per (Participant, Direction, ParameterSet,
%             Condition) POOLED group -- the finest level here (there is
%             no true single-trial level for this metric, by design).
%    Level B  Head_PCA_Deviation_GroupSummary.csv
%             Level A averaged ACROSS PARTICIPANTS per (ParameterSet,
%             Condition, Direction) -- exact-combo match, participant =
%             unit of replication.
%    Level C  Freq x Amp grid, same binning as the off-axis script:
%             amplitude via divisive tolerance clustering
%             (opts.amp_bin_tol_dps), frequency snapped to known nominal
%             designed values (opts.freq_nominal_targets):
%               Head_PCA_Deviation_FreqEffect_at_ConstAmplitude.csv
%               Head_PCA_Deviation_AmpEffect_at_ConstFrequency.csv
%             plus matching plots (x=frequency, lines=amplitude bin, and
%             vice versa; faceted by Direction, one figure set per
%             Condition present).
%    Level D  Head_PCA_Deviation_DirectionOverall_Summary.csv
%             One row per (Direction, Condition), collapsing across the
%             ENTIRE tested Freq x Amp space (equal weight per tested
%             input cell) -- the single "how much does the head's axis
%             deviate from the flange's, irrespective of input" number,
%             companion to the off-axis script's Level 4 (which uses
%             OffAxis_to_Task_RMS_Ratio as its headline magnitude
%             measure). Use both together: this angle is scale-invariant
%             (orientation only), the ratio is magnitude-based -- a
%             direction can have a small, stable angle but a ratio that
%             still grows with amplitude, or vice versa, and that
%             contrast is itself informative.
%
%  FOLDER LAYOUT ASSUMED (matches KUKA_analysis_preprocess_driver.m):
%     base_dir/P01_sync_summary_order/trimmed_*.mat
%     base_dir/P02_sync_summary_order/trimmed_*.mat
%     ...
%
%  For a per-Direction x per-Condition (EO/EC) BATCH -- one separate set
%  of these output files/plots per combination -- use the wrapper
%  KUKA_head_PCA_deviation_ByDirectionCondition.m.
%
%  USAGE
%    KUKA_head_PCA_deviation_AllParticipants(1:10)
%    filt = struct('Direction', 'ML', 'Condition', 'EO');
%    KUKA_head_PCA_deviation_AllParticipants(1:10, base_dir, filt)
% ========================================================================

%% ------------------------- CONFIG / DEFAULTS -----------------------------
if nargin < 1 || isempty(participant_nums)
    participant_nums = 1:10;
end
if nargin < 2 || isempty(base_dir)
    base_dir = pwd;
end
if nargin < 3 || isempty(filt)
    filt = struct();
end
filt = fill_default_filter(filt);
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = fill_default_opts(opts);

out_dir = opts.save_dir;
if isempty(out_dir)
    out_dir = fullfile(base_dir, 'AllParticipants_HeadPCADeviation');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ------------------------- PER-PARTICIPANT: DEDUP + POOL + PCA -----------
rows = {};

for pi = 1:numel(participant_nums)
    p = participant_nums(pi);
    csv_out_dir = participant_folder(base_dir, p);

    if ~exist(csv_out_dir, 'dir')
        fprintf('[SKIP participant %d] Folder not found: %s\n', p, csv_out_dir);
        continue;
    end
    mat_files = dir(fullfile(csv_out_dir, 'trimmed_*.mat'));
    if isempty(mat_files)
        fprintf('[SKIP participant %d] No trimmed_*.mat files in: %s\n', p, csv_out_dir);
        continue;
    end

    % --- Dedup: one file per PHYSICAL trial (lowest Muscle number), so
    % pooling below never double-counts a trial's cycles via its Muscle
    % 1/Muscle 2 duplicate .mat files (identical IMU data in both). ---
    best_by_trial = containers.Map();
    for k = 1:numel(mat_files)
        tok = parse_trimmed_filename(mat_files(k).name);
        if isempty(tok) || ~passes_filter(tok, filt)
            continue;
        end
        trial_key = strjoin({tok.imu_base, tok.direction, tok.paramset, tok.condition}, '|');
        entry = struct('name', mat_files(k).name, 'folder', mat_files(k).folder, 'tok', tok);
        if ~isKey(best_by_trial, trial_key) || tok.muscle < best_by_trial(trial_key).tok.muscle
            best_by_trial(trial_key) = entry;
        end
    end
    if best_by_trial.Count == 0
        fprintf('[SKIP participant %d] No files passed the filter.\n', p);
        continue;
    end
    dedup_list = values(best_by_trial);

    % --- Group the deduplicated physical trials by Direction+ParameterSet+Condition ---
    group_key_of = @(s) strjoin({s.tok.direction, s.tok.paramset, s.tok.condition}, '|');
    group_keys_all = cellfun(group_key_of, dedup_list, 'UniformOutput', false);
    [group_keys_unique, ~, gidx] = unique(group_keys_all);

    fprintf('\n=== Participant %d: %d physical trial(s), %d Direction/ParameterSet/Condition group(s) ===\n', ...
        p, numel(dedup_list), numel(group_keys_unique));

    for g = 1:numel(group_keys_unique)
        members = dedup_list(gidx == g);

        head_pool = [];
        flange_pool = [];
        files_included = {};
        for m = 1:numel(members)
            fpath = fullfile(members{m}.folder, members{m}.name);
            try
                S = load(fpath, 'data');
                data = S.data;
            catch ME
                fprintf('  [SKIP] Could not load %s: %s\n', members{m}.name, ME.message);
                continue;
            end
            if ~isfield(data, 'status') || ~strcmp(data.status, 'ok')
                continue;
            end
            if ~isfield(data, 'imu_head') || ~isfield(data.imu_head, 'six_axis_resamp') || ...
               ~isfield(data, 'imu_flange') || ~isfield(data.imu_flange, 'six_axis_resamp')
                continue;
            end
            head_pool   = [head_pool;   data.imu_head.six_axis_resamp(:, 4:6)];   %#ok<AGROW>
            flange_pool = [flange_pool; data.imu_flange.six_axis_resamp(:, 4:6)]; %#ok<AGROW>
            files_included{end+1} = members{m}.tok.imu_base; %#ok<AGROW>
        end

        if size(head_pool, 1) < 4 || size(flange_pool, 1) < 4
            fprintf('  [SKIP] %s: not enough pooled samples (%d).\n', group_keys_unique{g}, size(head_pool, 1));
            continue;
        end

        tok0 = members{1}.tok;  % Direction/ParameterSet/Condition are common across the group

        [head_PC1, head_VAF]     = pca_gyro(head_pool);
        [flange_PC1, flange_VAF] = pca_gyro(flange_pool);

        nominal = nominal_input_axis(tok0.direction);

        DeviationAngle_3D_deg       = acosd(min(1, abs(dot(head_PC1, flange_PC1))));
        Flange_OutOfPlane_Angle_deg = asind(min(1, abs(flange_PC1(3))));

        if any(isnan(nominal))
            head_PC1_sf = head_PC1;
        else
            head_PC1_sf = sign_fix_vec(head_PC1, nominal);
        end
        OutOfPlane_Angle_deg = asind(min(1, abs(head_PC1_sf(3))));

        % --- Per-CALIBRATED-AXIS (Gx,Gy,Gz) deviation vector -------------
        % head_PC1/flange_PC1 are already expressed in the calibrated
        % (Gx,Gy,Gz) frame from the two-pose calibration -- PCA doesn't
        % change coordinate systems, it just picks a direction WITHIN that
        % frame. head_PC1_sf above is sign-fixed against the NOMINAL design
        % axis (needed for OutOfPlane/InPlane above, which are defined
        % relative to the idealized direction), but a per-axis deviation
        % vector needs head sign-fixed against what the flange ACTUALLY
        % delivered this trial instead -- otherwise a diagonal trial whose
        % true achieved axis isn't exactly the nominal 45/-45 split would
        % have head and flange sign-referenced to two different things,
        % and the component-wise difference below would be meaningless.
        head_PC1_flangealigned = sign_fix_vec(head_PC1, flange_PC1);
        DeviationVec_Gx = head_PC1_flangealigned(1) - flange_PC1(1);
        DeviationVec_Gy = head_PC1_flangealigned(2) - flange_PC1(2);
        DeviationVec_Gz = head_PC1_flangealigned(3) - flange_PC1(3);

        head_inplane = [head_PC1_sf(1), head_PC1_sf(2)];
        n_ip = norm(head_inplane);
        if n_ip > 1e-6 && ~any(isnan(nominal))
            head_inplane_u = head_inplane / n_ip;
            nominal_ip = [nominal(1), nominal(2)];
            nominal_ip = nominal_ip / norm(nominal_ip);
            InPlane_Azimuth_Angle_deg = acosd(min(1, max(-1, dot(head_inplane_u, nominal_ip))));
            cross_z = nominal_ip(1) * head_inplane_u(2) - nominal_ip(2) * head_inplane_u(1);
            InPlane_Azimuth_SignedDeg = atan2d(cross_z, dot(nominal_ip, head_inplane_u));
        else
            InPlane_Azimuth_Angle_deg = NaN;
            InPlane_Azimuth_SignedDeg = NaN;
        end

        row = struct();
        row.Participant       = p;
        row.Direction          = string(tok0.direction);
        row.ParameterSet       = string(tok0.paramset);
        row.Freq_Hz            = tok0.freq_hz;
        row.Amp_dps             = tok0.amp_dps;
        row.Condition           = string(tok0.condition);
        row.N_Trials_Pooled      = numel(files_included);
        row.N_Samples_Pooled      = size(head_pool, 1);
        row.TrialFiles_Included    = strjoin(files_included, ',');
        row.Head_PC1_VAF            = head_VAF;
        row.Flange_PC1_VAF           = flange_VAF;
        row.Flange_OutOfPlane_Angle_deg = Flange_OutOfPlane_Angle_deg;
        row.DeviationAngle_3D_deg        = DeviationAngle_3D_deg;
        row.OutOfPlane_Angle_deg          = OutOfPlane_Angle_deg;
        row.InPlane_Azimuth_Angle_deg      = InPlane_Azimuth_Angle_deg;
        row.InPlane_Azimuth_SignedDeg       = InPlane_Azimuth_SignedDeg;
        row.Head_PC1_x = head_PC1_sf(1); row.Head_PC1_y = head_PC1_sf(2); row.Head_PC1_z = head_PC1_sf(3);
        row.Flange_PC1_x = flange_PC1(1); row.Flange_PC1_y = flange_PC1(2); row.Flange_PC1_z = flange_PC1(3);
        row.DeviationVec_Gx = DeviationVec_Gx;
        row.DeviationVec_Gy = DeviationVec_Gy;
        row.DeviationVec_Gz = DeviationVec_Gz;

        rows{end+1} = row; %#ok<AGROW>

        fprintf('  %-30s : %d trial(s), %d samples pooled -- DeviationAngle=%.1f deg (Head VAF=%.2f, Flange VAF=%.2f)\n', ...
            group_keys_unique{g}, numel(files_included), size(head_pool,1), DeviationAngle_3D_deg, head_VAF, flange_VAF);
    end
end

if isempty(rows)
    error(['No usable (Participant, Direction, ParameterSet, Condition) groups found under:\n  %s\n' ...
        'Check base_dir, participant_nums, and filt.'], base_dir);
end

pooled_table = struct2table([rows{:}]);

%% ------------------------- LEVEL A: PER-PARTICIPANT POOLED TABLE ---------
out_csv_A = fullfile(out_dir, 'Head_PCA_Deviation_PerParticipant.csv');
writetable(pooled_table, out_csv_A);
fprintf('\nLevel A (per-participant, pooled-trial PCA) written to:\n  %s\n', out_csv_A);

%% ------------------------- LEVEL B: GROUP SUMMARY (exact combo) ----------
lb_keys = strcat(pooled_table.ParameterSet, "_", pooled_table.Condition, "_", pooled_table.Direction);
[lb_unique, ~, lb_idx] = unique(lb_keys);

is_num = varfun(@isnumeric, pooled_table, 'OutputFormat', 'uniform');
numeric_vars = pooled_table.Properties.VariableNames(is_num);
numeric_vars = setdiff(numeric_vars, {'Participant'}, 'stable');

lb_rows = {};
for g = 1:numel(lb_unique)
    idx = (lb_idx == g);
    sub = pooled_table(idx, :);
    srow = struct();
    srow.ParameterSet   = sub.ParameterSet(1);
    srow.Condition      = sub.Condition(1);
    srow.Direction      = sub.Direction(1);
    srow.Freq_Hz        = sub.Freq_Hz(1);
    srow.Amp_dps        = sub.Amp_dps(1);
    srow.N_Participants = height(sub);
    srow.Participants_Included = strjoin(string(unique(sub.Participant)), ',');
    for v = 1:numel(numeric_vars)
        vn = numeric_vars{v};
        vals = sub.(vn);
        m = mean(vals, 'omitnan');
        s = std(vals, 'omitnan');
        srow.([vn '_mean']) = m;
        srow.([vn '_between_SD']) = s;
        if m ~= 0
            srow.([vn '_between_CV']) = s / m;
        else
            srow.([vn '_between_CV']) = NaN;
        end
    end
    lb_rows{end+1} = srow; %#ok<AGROW>
end
lb_table = struct2table([lb_rows{:}]);
lb_table = sortrows(lb_table, {'Direction', 'Freq_Hz', 'Amp_dps', 'Condition'});

out_csv_B = fullfile(out_dir, 'Head_PCA_Deviation_GroupSummary.csv');
writetable(lb_table, out_csv_B);
fprintf('Level B (group summary across participants, exact combo) written to:\n  %s\n', out_csv_B);
disp(lb_table(:, {'Direction','Freq_Hz','Amp_dps','Condition','N_Participants', ...
    'DeviationAngle_3D_deg_mean','DeviationAngle_3D_deg_between_SD'}));

%% ------------------------- LEVEL C: FREQ x AMP GRID -----------------------
% AMPLITUDE: divisive gap-based tolerance clustering (opts.amp_bin_tol_dps).
% FREQUENCY: snapped to the known nominal designed values
% (opts.freq_nominal_targets) instead -- see the off-axis script's Level C
% section for the full rationale. Set opts.freq_nominal_targets = [] to
% fall back to tolerance clustering via opts.freq_bin_tol_hz.
amp_tol = opts.amp_bin_tol_dps;
[amp_bin_id, amp_bin_centers] = cluster_values_tolerance(pooled_table.Amp_dps, amp_tol);
pooled_table.AmpBin = round(amp_bin_centers(amp_bin_id), 3);

if isempty(opts.freq_nominal_targets)
    freq_tol = opts.freq_bin_tol_hz;
    [freq_bin_id, freq_bin_centers] = cluster_values_tolerance(pooled_table.Freq_Hz, freq_tol);
    pooled_table.FreqBin = round(freq_bin_centers(freq_bin_id), 4);
else
    nominal_targets = opts.freq_nominal_targets(:);
    freq_bin_id = snap_to_nominal(pooled_table.Freq_Hz, nominal_targets);
    pooled_table.FreqBin = nominal_targets(freq_bin_id);
end

grid_keys = strcat(pooled_table.Direction, "_", pooled_table.Condition, "_", ...
    string(pooled_table.FreqBin), "_", string(pooled_table.AmpBin));
[grid_unique, ~, grid_idx] = unique(grid_keys);

is_num2 = varfun(@isnumeric, pooled_table, 'OutputFormat', 'uniform');
numeric_vars2 = pooled_table.Properties.VariableNames(is_num2);
numeric_vars2 = setdiff(numeric_vars2, {'Participant', 'FreqBin', 'AmpBin'}, 'stable');

grid_rows = {};
for g = 1:numel(grid_unique)
    idx = (grid_idx == g);
    sub = pooled_table(idx, :);
    srow = struct();
    srow.Direction   = sub.Direction(1);
    srow.Condition    = sub.Condition(1);
    srow.FreqBin       = sub.FreqBin(1);
    srow.AmpBin          = sub.AmpBin(1);
    srow.N_Participants    = height(sub);
    srow.Participants_Included = strjoin(string(unique(sub.Participant)), ',');
    srow.ParameterSets_Included = strjoin(unique(sub.ParameterSet), ',');
    for v = 1:numel(numeric_vars2)
        vn = numeric_vars2{v};
        vals = sub.(vn);
        m = mean(vals, 'omitnan');
        s = std(vals, 'omitnan');
        srow.([vn '_mean']) = m;
        srow.([vn '_between_SD']) = s;
        if m ~= 0
            srow.([vn '_between_CV']) = s / m;
        else
            srow.([vn '_between_CV']) = NaN;
        end
    end
    grid_rows{end+1} = srow; %#ok<AGROW>
end
grid_table = struct2table([grid_rows{:}]);

out_csv_freq = fullfile(out_dir, 'Head_PCA_Deviation_FreqEffect_at_ConstAmplitude.csv');
writetable(sortrows(grid_table, {'Direction', 'Condition', 'AmpBin', 'FreqBin'}), out_csv_freq);
fprintf('Level C frequency-effect table written to:\n  %s\n', out_csv_freq);

out_csv_amp = fullfile(out_dir, 'Head_PCA_Deviation_AmpEffect_at_ConstFrequency.csv');
writetable(sortrows(grid_table, {'Direction', 'Condition', 'FreqBin', 'AmpBin'}), out_csv_amp);
fprintf('Level C amplitude-effect table written to:\n  %s\n', out_csv_amp);

%% ------------------- SHARED-ACROSS-DIRECTIONS FLAG -----------------------
% See the off-axis script's copy of this block for the full rationale:
% not every direction could achieve the same designed frequency/
% amplitude, so a cell is only "fair game" for cross-direction comparison
% if every tested direction actually has data in it.
grid_table.SharedAcrossDirections = false(height(grid_table), 1);
conditions_present = unique(grid_table.Condition, 'stable');
for ci = 1:numel(conditions_present)
    cond = conditions_present(ci);
    cond_mask = (grid_table.Condition == cond);
    n_dir_total = numel(unique(grid_table.Direction(cond_mask)));

    cell_keys_all = strcat(string(grid_table.FreqBin), "_", string(grid_table.AmpBin));
    rows_in_cond = find(cond_mask);
    [cell_unique, ~, cell_idx] = unique(cell_keys_all(rows_in_cond));
    for ck = 1:numel(cell_unique)
        members = rows_in_cond(cell_idx == ck);
        n_dir_here = numel(unique(grid_table.Direction(members)));
        if n_dir_here == n_dir_total
            grid_table.SharedAcrossDirections(members) = true;
        end
    end
end

%% ------------------------- LEVEL D: DIRECTION/CONDITION OVERALL ----------
dc_keys = strcat(grid_table.Direction, "_", grid_table.Condition);
[dc_unique, ~, dc_idx] = unique(dc_keys);

metric_cols = { ...
    'DeviationAngle_3D_deg_mean',      '3D head-vs-flange PCA deviation angle (deg) -- HEADLINE'; ...
    'OutOfPlane_Angle_deg_mean',       'Out-of-plane (yaw) component of head axis (deg)'; ...
    'InPlane_Azimuth_Angle_deg_mean',  'In-plane azimuthal bearing error (deg)'; ...
    'DeviationVec_Gx_mean',            'Per-calibrated-axis deviation, Gx component (unitless, flange-aligned)'; ...
    'DeviationVec_Gy_mean',            'Per-calibrated-axis deviation, Gy component (unitless, flange-aligned)'; ...
    'DeviationVec_Gz_mean',            'Per-calibrated-axis deviation, Gz component (unitless, flange-aligned)'; ...
    'Head_PC1_VAF_mean',               'Head PC1 variance-explained (reliability of the angle above)'};

d4_rows = {};
for g = 1:numel(dc_unique)
    idx = (dc_idx == g);
    sub = grid_table(idx, :);
    drow = struct();
    drow.Direction    = sub.Direction(1);
    drow.Condition     = sub.Condition(1);
    drow.N_FreqAmpCells  = height(sub);
    match = (pooled_table.Direction == sub.Direction(1)) & (pooled_table.Condition == sub.Condition(1));
    drow.N_Participants_Total = numel(unique(pooled_table.Participant(match)));
    for mi = 1:size(metric_cols, 1)
        col = metric_cols{mi, 1};
        if ~ismember(col, sub.Properties.VariableNames)
            continue;
        end
        vals = sub.(col);
        m = mean(vals, 'omitnan');
        s = std(vals, 'omitnan');
        drow.([col '_GrandMean']) = m;
        drow.([col '_AcrossInputSD']) = s;
        if m ~= 0
            drow.([col '_AcrossInputCV']) = s / m;
        else
            drow.([col '_AcrossInputCV']) = NaN;
        end
    end
    d4_rows{end+1} = drow; %#ok<AGROW>
end
direction_summary_table = struct2table([d4_rows{:}]);
direction_summary_table = sortrows(direction_summary_table, {'Direction', 'Condition'});

out_csv_D = fullfile(out_dir, 'Head_PCA_Deviation_DirectionOverall_Summary.csv');
writetable(direction_summary_table, out_csv_D);
fprintf(['Level D (direction/condition overall summary, irrespective of input) written to:\n  %s\n' ...
    '  NOTE: not directly comparable direction-to-direction -- see\n' ...
    '  Head_PCA_Deviation_DirectionComparison_SharedBinsOnly.csv below.\n'], out_csv_D);
disp(direction_summary_table);

%% --------- LEVEL D-SHARED: DIRECTION COMPARISON (SHARED BINS ONLY) --------
shared_grid = grid_table(grid_table.SharedAcrossDirections, :);
if isempty(shared_grid)
    fprintf('\n[WARNING] No (FreqBin, AmpBin) cell is shared across every tested direction -- Level D-Shared skipped.\n');
else
    dcs_keys = strcat(shared_grid.Direction, "_", shared_grid.Condition);
    [dcs_unique, ~, dcs_idx] = unique(dcs_keys);

    d5_rows = {};
    for g = 1:numel(dcs_unique)
        idx = (dcs_idx == g);
        sub = shared_grid(idx, :);

        drow = struct();
        drow.Direction              = sub.Direction(1);
        drow.Condition               = sub.Condition(1);
        drow.N_SharedFreqAmpCells    = height(sub);
        drow.SharedFreqBins_Included = strjoin(string(unique(sub.FreqBin)), ',');
        drow.SharedAmpBins_Included  = strjoin(string(unique(sub.AmpBin)), ',');

        pooled_cell_keys = strcat(string(pooled_table.FreqBin), "_", string(pooled_table.AmpBin));
        sub_cell_keys    = strcat(string(sub.FreqBin), "_", string(sub.AmpBin));
        match = (pooled_table.Direction == sub.Direction(1)) & (pooled_table.Condition == sub.Condition(1)) & ...
            ismember(pooled_cell_keys, sub_cell_keys);
        drow.N_Participants_Total = numel(unique(pooled_table.Participant(match)));

        for mi = 1:size(metric_cols, 1)
            col = metric_cols{mi, 1};
            if ~ismember(col, sub.Properties.VariableNames)
                continue;
            end
            vals = sub.(col);
            m = mean(vals, 'omitnan');
            s = std(vals, 'omitnan');
            drow.([col '_GrandMean']) = m;
            drow.([col '_AcrossInputSD']) = s;
            if m ~= 0
                drow.([col '_AcrossInputCV']) = s / m;
            else
                drow.([col '_AcrossInputCV']) = NaN;
            end
        end
        d5_rows{end+1} = drow; %#ok<AGROW>
    end
    direction_comparison_table = struct2table([d5_rows{:}]);
    direction_comparison_table = sortrows(direction_comparison_table, {'Condition', 'Direction'});

    out_csv_Dshared = fullfile(out_dir, 'Head_PCA_Deviation_DirectionComparison_SharedBinsOnly.csv');
    writetable(direction_comparison_table, out_csv_Dshared);
    fprintf('Level D-Shared (fair cross-direction comparison, shared Freq/Amp cells only) written to:\n  %s\n', out_csv_Dshared);
    disp(direction_comparison_table);
end

%% ------------- DIAGNOSTIC: (FreqBin, AmpBin) CELLS NOT SHARED --------------
excluded_grid = grid_table(~grid_table.SharedAcrossDirections, :);
if ~isempty(excluded_grid)
    out_csv_excluded = fullfile(out_dir, 'Head_PCA_Deviation_NonSharedFreqAmpCells.csv');
    writetable(sortrows(excluded_grid(:, {'Direction', 'Condition', 'FreqBin', 'AmpBin', 'N_Participants'}), ...
        {'Condition', 'FreqBin', 'AmpBin', 'Direction'}), out_csv_excluded);
    fprintf(['\n%d (FreqBin, AmpBin) cell(s) are NOT shared across every tested direction\n' ...
        '(excluded from Level D-Shared) -- see:\n  %s\n'], height(excluded_grid), out_csv_excluded);
end

%% ------------------------- PLOTS ------------------------------------------
plot_pca_deviation_effect(grid_table, out_dir, 'FreqBin', 'AmpBin', ...
    'Input frequency (Hz)', 'Amp %.2f dps', 'PCADeviation_FreqEffect_at_ConstAmplitude');
plot_pca_deviation_effect(grid_table, out_dir, 'AmpBin', 'FreqBin', ...
    'Input amplitude (deg/s)', 'Freq %.3f Hz', 'PCADeviation_AmpEffect_at_ConstFrequency');

end % KUKA_head_PCA_deviation_AllParticipants


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function csv_out_dir = participant_folder(base_dir, p)
    csv_out_dir = fullfile(base_dir, sprintf('P%02d_sync_summary_order', p));
end

function filt = fill_default_filter(filt)
    if ~isfield(filt, 'ParameterSet'); filt.ParameterSet = ''; end
    if ~isfield(filt, 'Condition');    filt.Condition    = ''; end
    if ~isfield(filt, 'Direction');    filt.Direction    = ''; end
end

function opts = fill_default_opts(opts)
    if ~isfield(opts, 'save_dir');           opts.save_dir = '';       end
    if ~isfield(opts, 'freq_bin_tol_hz');    opts.freq_bin_tol_hz = 0.02; end
    if ~isfield(opts, 'amp_bin_tol_dps');    opts.amp_bin_tol_dps = 0.2;  end
    if ~isfield(opts, 'freq_nominal_targets')
        % Known DESIGNED frequency levels (Hz) -- see the off-axis
        % script's fill_default_opts for the full rationale. Set to []
        % to fall back to tolerance clustering via freq_bin_tol_hz.
        opts.freq_nominal_targets = [0.4, 0.48, 0.55, 0.8, 1.0, 1.3];
        % 0.48 kept separate from 0.55 -- see the off-axis script's copy
        % of this default for the full rationale (diagonal-only
        % compromise value, distinct from ML's true 0.55 Hz).
    end
end

function bin_id = snap_to_nominal(vals, nominal_targets)
% Assigns each value in vals to whichever entry in nominal_targets it is
% numerically closest to. See the off-axis script's copy of this
% function for the full explanation.
    vals = vals(:);
    nominal_targets = nominal_targets(:);
    diffs = abs(vals - nominal_targets');
    [~, bin_id] = min(diffs, [], 2);
end

function tok = parse_trimmed_filename(fname)
% ParameterSet = "<freq*1000>_<amp*100>", e.g. "1800_0403" ->
%   freq_hz = 1.800 Hz, amp_raw = 403, amp_dps = 4.03 deg/s
    tok = [];
    pat = '^trimmed_(.+)_Muscle(\d+)_([A-Za-z]+)_(\d+_\d+)_([A-Za-z0-9]+)\.mat$';
    m = regexp(fname, pat, 'tokens', 'once');
    if isempty(m)
        return;
    end
    tok = struct();
    tok.imu_base  = m{1};
    tok.muscle    = str2double(m{2});
    tok.direction = m{3};
    tok.paramset  = m{4};
    tok.condition = m{5};

    freq_tok = regexp(tok.paramset, '^(\d+)_(\d+)$', 'tokens', 'once');
    if ~isempty(freq_tok)
        tok.freq_hz = str2double(freq_tok{1}) / 1000;
        tok.amp_raw = str2double(freq_tok{2});
        tok.amp_dps = tok.amp_raw / 100;
    else
        tok.freq_hz = NaN;
        tok.amp_raw = NaN;
        tok.amp_dps = NaN;
    end
end

function tf = passes_filter(tok, filt)
    tf = true;
    if ~isempty(filt.ParameterSet) && ~strcmpi(tok.paramset, filt.ParameterSet)
        tf = false; return;
    end
    if ~isempty(filt.Condition) && ~strcmpi(tok.condition, filt.Condition)
        tf = false; return;
    end
    if ~isempty(filt.Direction) && ~strcmpi(tok.direction, filt.Direction)
        tf = false; return;
    end
end

function [bin_id, bin_center] = cluster_values_tolerance(vals, tol)
% DIVISIVE gap-based clustering (same as the off-axis script): looks at
% the whole sorted set at once, starts with everything in one bin, and
% repeatedly splits any bin whose width (max-min) exceeds tol at its OWN
% LARGEST INTERNAL GAP, until every bin is <= tol wide. This guarantees
% no bin ever exceeds tol (no unbounded chaining) while avoiding the
% earlier from-bin-start version's issue of splitting two values that
% are close to EACH OTHER (e.g. 4.02/4.09 at tol=0.2) just because the
% bin happened to be anchored at some unrelated lower value. See the
% off-axis script's copy of this function for the fuller explanation.
    vals = vals(:);
    n = numel(vals);
    [sorted_vals, sort_idx] = sort(vals);

    tol_eps = 1e-9;   % guards against float rounding spuriously splitting a
                      % group whose true span sits exactly at tol
    clusters = {[1, n]};
    changed = true;
    while changed
        changed = false;
        next_clusters = {};
        for c = 1:numel(clusters)
            lo = clusters{c}(1);
            hi = clusters{c}(2);
            if hi > lo && (sorted_vals(hi) - sorted_vals(lo)) > tol + tol_eps
                seg_gaps = diff(sorted_vals(lo:hi));
                [~, gi] = max(seg_gaps);
                split_at = lo + gi - 1;
                next_clusters{end+1} = [lo, split_at];       %#ok<AGROW>
                next_clusters{end+1} = [split_at + 1, hi];   %#ok<AGROW>
                changed = true;
            else
                next_clusters{end+1} = clusters{c}; %#ok<AGROW>
            end
        end
        clusters = next_clusters;
    end

    starts = cellfun(@(c) c(1), clusters);
    [~, order] = sort(starts);
    clusters = clusters(order);

    bin_id_sorted = nan(n, 1);
    for b = 1:numel(clusters)
        lo = clusters{b}(1);
        hi = clusters{b}(2);
        bin_id_sorted(lo:hi) = b;
    end

    bin_id = nan(n, 1);
    bin_id(sort_idx) = bin_id_sorted;

    n_bins = numel(clusters);
    bin_center = nan(n_bins, 1);
    for b = 1:n_bins
        bin_center(b) = mean(vals(bin_id == b), 'omitnan');
    end
end

function [PC1, VAF] = pca_gyro(gyro_data)
% 3D PCA on a pooled [N x 3] (Gx,Gy,Gz) sample: mean-centers, computes the
% 3x3 covariance, and returns the dominant eigenvector (unit vector) and
% the fraction of total variance it explains (VAF -- a reliability
% indicator: near 1 = well-defined single axis, near 1/3 = isotropic/
% undefined axis).
    gyro_data = gyro_data - mean(gyro_data, 1, 'omitnan');
    C = cov(gyro_data);
    [V, D] = eig(C);
    eigvals = diag(D);
    [eigvals_sorted, order] = sort(eigvals, 'descend');
    V = V(:, order);
    PC1 = V(:, 1);
    PC1 = PC1 / norm(PC1);
    total_var = sum(eigvals_sorted);
    if total_var > 0
        VAF = eigvals_sorted(1) / total_var;
    else
        VAF = NaN;
    end
end

function v = nominal_input_axis(direction)
% Same nominal per-direction reference convention already used elsewhere
% in this pipeline (e.g. the diagonal-trial PCA sign-fix in
% KUKA_analysis_prep_consume_ref.m Section 5 / reconstruct_flange_reference).
    d = upper(direction);
    if strcmp(d, 'ML')
        v = [1; 0; 0];
    elseif strcmp(d, 'AP')
        v = [0; 1; 0];
    elseif contains(d, 'DIA')
        if contains(d, 'R')
            v = [1; 1; 0];
        else
            v = [1; -1; 0];
        end
        v = v / norm(v);
    else
        v = [NaN; NaN; NaN];  % unrecognized direction -- don't assume an axis
    end
end

function v = sign_fix_vec(v, ref)
% Flips v (a unit-norm eigenvector, only defined up to +/-1) so it points
% into the same half-space as ref.
    if dot(v, ref) < 0
        v = -v;
    end
end

function plot_pca_deviation_effect(grid_table, out_dir, x_field, line_field, x_label, line_label_fmt, file_prefix)
% Same structure as the off-axis script's plot_offaxis_effect, but for
% the PCA deviation metrics, and split into one figure set PER CONDITION
% present (since grid_table's Direction subplots would otherwise mix
% EO/EC points at the same FreqBin/AmpBin together).
    try
        conditions = unique(grid_table.Condition, 'stable');
        directions = unique(grid_table.Direction, 'stable');
        metrics = { ...
            'DeviationAngle_3D_deg',   '3D deviation angle (deg)',    '_Deviation3D.png'; ...
            'OutOfPlane_Angle_deg',    'Out-of-plane (yaw) angle (deg)', '_OutOfPlane.png' };

        for ci = 1:numel(conditions)
            cond = conditions(ci);
            csub = grid_table(grid_table.Condition == cond, :);
            if isempty(csub)
                continue;
            end

            for mi = 1:size(metrics, 1)
                mean_col = [metrics{mi,1} '_mean'];
                sd_col   = [metrics{mi,1} '_between_SD'];
                if ~ismember(mean_col, csub.Properties.VariableNames)
                    continue;
                end

                fig = figure('Visible', 'off', 'Position', [100 100 1400 350 * numel(directions)]);
                for di = 1:numel(directions)
                    subplot(numel(directions), 1, di);
                    dsub = csub(csub.Direction == directions(di), :);
                    line_vals = unique(dsub.(line_field));
                    hold on;
                    for li = 1:numel(line_vals)
                        lsub = dsub(dsub.(line_field) == line_vals(li), :);
                        lsub = sortrows(lsub, x_field);
                        errorbar(lsub.(x_field), lsub.(mean_col), lsub.(sd_col), '-o', ...
                            'DisplayName', sprintf(line_label_fmt, line_vals(li)));
                    end
                    hold off;
                    grid on;
                    xlabel(x_label);
                    ylabel(metrics{mi,2});
                    title(sprintf('Direction: %s | Condition: %s', directions(di), cond));
                    legend('Location', 'best');
                end
                sgtitle(strrep(metrics{mi,2}, '_', '\_'));
                out_png = fullfile(out_dir, [file_prefix '_' char(cond) metrics{mi,3}]);
                print(fig, out_png, '-dpng', '-r150');
                close(fig);
                fprintf('Saved: %s\n', out_png);
            end
        end
    catch ME
        fprintf('[WARNING] Could not generate %s plots: %s\n', file_prefix, ME.message);
    end
end
