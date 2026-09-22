function KUKA_head_kinematics_summary_AllParticipants(participant_nums, base_dir, filt, opts)
%% ========================================================================
%  HEAD KINEMATICS + OFF-AXIS COMPENSATION ERROR -- ALL-PARTICIPANT BATCH
% ========================================================================
%  Multi-participant wrapper around the per-trial kinematics logic in
%  KUKA_head_kinematics_summary_v4.m. Loops every participant's
%  trimmed_*.mat files (written by KUKA_analysis_preprocess_consumer.m,
%  driven by KUKA_analysis_preprocess_driver.m), PLUS a magnitude+
%  variability off-axis deviation bundle (see "OFF-AXIS DEVIATION METRICS"
%  below), then aggregates across FIVE levels so frequency/amplitude
%  effects can be examined without pseudoreplication.
%
%  *** TRIALS SHARING THE EXACT SAME CONDITION ARE POOLED BEFORE ANY
%  METRIC IS COMPUTED, NOT AVERAGED AFTER. *** One (Direction, Freq, Amp,
%  Condition) input condition is typically delivered as 40 total cycles,
%  which may be split across trials differently (e.g. 2 trials of 20
%  cycles, or 4 trials of 10) -- if you computed RMS/PtP/TotalExcursion/
%  etc. on each trial separately and then just averaged those numbers,
%  the result would implicitly depend on how the 40 cycles happened to be
%  split into trials (a 4-trial condition would weight each 10-cycle
%  chunk 1/4, a 2-trial condition would weight each 20-cycle chunk 1/2,
%  even though both represent the same total data), and quantities like
%  TotalExcursion (a cumulative trapz over time) aren't even comparable
%  across trials of different duration in the first place. So instead:
%  every physical trial sharing the same (Participant, Direction,
%  ParameterSet, Condition) has its cycle-trimmed, resampled six-axis
%  gyro/accel samples (data.imu_head.six_axis_resamp / data.imu_flange.
%  six_axis_resamp) CONCATENATED
%  end-to-end first, and every metric below (RMS, PtP, Bias,
%  TotalExcursion, PhaseCorr/Lag, the off-axis resultant bundle, ...) is
%  computed ONCE on that pooled sample. This naturally weights each trial
%  by how many samples/cycles it actually contributed, independent of how
%  the total was split across trials.
%
%  Each physical trial is also written out once per Muscle bank (1 and 2)
%  with IDENTICAL head/flange IMU data (only the EMG channels differ) --
%  those duplicates are collapsed to one file per physical trial (lowest
%  Muscle number kept) BEFORE pooling, so a trial's cycles are never
%  double-counted just because both its Muscle-bank files exist.
%
%  OUTPUT LEVELS
%    Level 1 (diagnostic only, NOT used by anything downstream)
%             Head_Kinematics_PerTrial_Diagnostic.csv
%             One row per individual physical trial, computed the OLD way
%             (metrics on that one trial alone) purely so you can eyeball
%             an individual trial (e.g. to spot an outlier) before
%             trusting the pooled result. Every aggregate below is built
%             from the POOLED computation, not from averaging these rows.
%
%    Level A  Head_Kinematics_PerParticipant_Pooled.csv
%             One row per (Participant, Direction, ParameterSet,
%             Condition) -- the finest level actually used for
%             aggregation. Metrics computed ONCE on the pooled/
%             concatenated samples from every (deduplicated) physical
%             trial sharing that exact condition, for that participant.
%
%    Level B  Head_Kinematics_Group_Summary.csv   <-- exact-combo freq/amp table
%             Level A averaged ACROSS PARTICIPANTS per (ParameterSet,
%             Condition, Direction) -- exact-combo match, participant is
%             the unit of replication. Reports both the group MEAN and
%             the BETWEEN-PARTICIPANT SD/CV of every numeric column.
%
%    Level C  FREQ x AMP GRID -- for the independent frequency/amplitude
%             effect question: same mean/between-SD/CV logic as Level B,
%             but nearby-but-not-identical designed combinations are
%             first collapsed into one FreqBin / AmpBin before grouping.
%             AMPLITUDE uses divisive gap-based clustering (Amp_dps
%             values within +/- opts.amp_bin_tol_dps, default 0.2 deg/s,
%             of each other collapse together -- validated against this
%             dataset's actual amplitude gaps). FREQUENCY is instead
%             SNAPPED TO THE KNOWN NOMINAL DESIGNED VALUES in
%             opts.freq_nominal_targets (default [0.4 0.55 0.8 1.0 1.3]
%             Hz) -- every observed Freq_Hz is assigned to whichever
%             target it's numerically closest to, since frequency's real
%             jitter doesn't have one clean universal gap size the way
%             amplitude does (set opts.freq_nominal_targets = [] to fall
%             back to tolerance clustering via opts.freq_bin_tol_hz
%             instead, e.g. if your design's frequencies change). Written
%             out as TWO CSVs (same underlying grid, reorganized for each
%             question) plus two pairs of plots:
%
%               Head_Kinematics_FreqEffect_at_ConstAmplitude.csv
%                 sorted Direction, Condition, AmpBin, FreqBin -- i.e.
%                 "holding amplitude ~constant, how do magnitude/
%                 variance change across frequency."
%               OffAxis_FreqEffect_at_ConstAmplitude_Magnitude_<Condition>.png
%               OffAxis_FreqEffect_at_ConstAmplitude_Variability_<Condition>.png
%                 (one PNG PER CONDITION present in the data -- e.g. _EO
%                 and _EC written separately, never overlaid on the same
%                 axes/line, since EO and EC are a comparison of interest
%                 and mixing them onto one line would be misleading)
%
%               Head_Kinematics_AmpEffect_at_ConstFrequency.csv
%                 sorted Direction, Condition, FreqBin, AmpBin -- i.e.
%                 "holding frequency ~constant, how do magnitude/
%                 variance change across amplitude."
%               OffAxis_AmpEffect_at_ConstFrequency_Magnitude_<Condition>.png
%               OffAxis_AmpEffect_at_ConstFrequency_Variability_<Condition>.png
%                 (also one PNG per condition, same reasoning as above)
%
%    Level D  Head_Kinematics_DirectionOverall_Summary.csv
%             One row per (Direction, Condition) present in this run,
%             collapsing across the ENTIRE tested Freq x Amp space for
%             THAT direction -- "how much off-axis movement happens in
%             this direction, irrespective of input." Headline metric is
%             OffAxis_to_Task_RMS_Ratio_Resamp_mean (off-axis RMS
%             normalized by task-axis RMS -- NOT raw off-axis RMS, which
%             is confounded with how hard the input itself was driving
%             movement). CAUTION: not all directions could achieve the
%             same designed frequency/amplitude (a compromise value was
%             substituted for some), so different directions' rows here
%             are integrated over DIFFERENT input spaces -- this table is
%             NOT meant for direction-to-direction comparison. Use the
%             next output for that.
%
%    Level D-Shared  Head_Kinematics_DirectionComparison_SharedBinsOnly.csv
%             Same grand-mean/AcrossInput-SD/CV calculation as Level D,
%             but restricted to only the (FreqBin, AmpBin) grid cells
%             that have data for EVERY direction tested in that condition
%             (grid_table.SharedAcrossDirections) -- a fair like-for-like
%             basis for comparing directions against each other. Cells
%             that exist for only some directions (e.g. AP's ~3.6/~6.9
%             deg/s amplitudes, or the diagonals' ~0.48 Hz frequency, that
%             no other direction ever tested) are left out of this table
%             and listed instead in
%             Head_Kinematics_NonSharedFreqAmpCells.csv, so you can see
%             exactly what was excluded and why.
%
%  For a per-Direction x per-Condition (EO/EC) BATCH -- i.e. one separate
%  set of these output files/plots per Direction x Condition combination
%  instead of one combined run -- use the wrapper
%  KUKA_head_kinematics_summary_ByDirectionCondition.m, which just calls
%  this function once per combination with filt/opts.save_dir set
%  appropriately.
%
%  FOLDER LAYOUT ASSUMED (matches KUKA_analysis_preprocess_driver.m):
%     base_dir/P01_sync_summary_order/trimmed_*.mat
%     base_dir/P02_sync_summary_order/trimmed_*.mat
%     ...
%
%  AMPLITUDE UNITS: the ParameterSet string's trailing number is in
%  centi-deg/s -- e.g. "0403" in "1800_0403" means 4.03 deg/s -- so
%  Amp_dps = Amp_raw / 100 throughout this file.
%
%  OFF-AXIS DEVIATION METRICS
%  ---------------------------------------
%  For every axis (Gx, Gy, Gz), this script computes RMS/PtP/
%  TotalExcursion/Bias/PhaseCorr/PhaseLag on the pooled sample, plus a
%  RESULTANT off-axis bundle: the off-axis gyro channel(s) are combined
%  into one instantaneous magnitude time series --
%       OffAxisResultant(t) = sqrt( sum over off-axis gyro channels of
%                                    signal(t)^2 )
%  (e.g. sqrt(Gy(t)^2 + Gz(t)^2) for ML; just |Gz(t)| for a diagonal
%  trial) -- respecting that the off-axis channels are simultaneous/
%  correlated, not independent noise sources. From that pooled time
%  series:
%    OffAxis_Resultant_RMS/Mean/PtP/TotalExcursion   magnitude
%    OffAxis_Resultant_SD / CV (=SD/Mean)             variability
%
%  At Level B / Level C, every _mean column above also gets a paired
%  _between_SD / _between_CV column computed ACROSS PARTICIPANTS.
%
%  PER-AXIS GYRO ENERGY FRACTION (Head_GyroEnergyFrac_Gx/Gy/Gz)
%  ---------------------------------------
%  Head_Gx/Gy/Gz_RMS already give the per-axis magnitude individually, but
%  OffAxis_Resultant_RMS/CV above combine whichever axes are off-axis for
%  a given direction into one number, which can obscure WHICH off-axis
%  channel actually dominates. Head_GyroEnergyFrac_Gx/Gy/Gz report what
%  FRACTION of the head's total rotational energy (sum of squared RMS
%  across all 3 axes) landed on each calibrated axis -- unsigned,
%  magnitude-based (consistent with RMS elsewhere in this script), always
%  summing to 1, and defined the same way regardless of direction (no
%  per-direction special-casing of which axes count as "off-axis"). E.g.
%  an ML trial might show ~0.80 Gx (task), 0.15 Gy, 0.05 Gz. This is the
%  magnitude-domain counterpart to KUKA_head_PCA_deviation_
%  AllParticipants.m's DeviationVec_Gx/Gy/Gz, which is signed/directional
%  (head axis vs. flange axis) rather than an energy fraction.
%
%  WITHIN-TRIAL vs BETWEEN-TRIAL DECOMPOSITION (Level A, per participant)
%  ------------------------------------------------------------------
%  OffAxis_Resultant_SD/_Mean (above) come from concatenating every trial
%  in a condition set and computing one variance/mean -- this is the
%  mathematically exact TOTAL variance/mean of the condition set, but it
%  is a blend of two different sources: how much the signal varies
%  cycle-to-cycle WITHIN a trial, and how much trials' own means differ
%  from each other (drift, block effects, or genuine behavioral
%  difference between exposures to the same nominal condition). To
%  separate these (law of total variance, one-way random-effects form),
%  each pooled Resamp block also reports, per trial i with n_i
%  samples, mean m_i, variance v_i, and weighted grand mean
%  gm = sum(n_i*m_i)/N:
%    OffAxis_Resultant_WithinTrial_SD        = sqrt( sum(n_i*v_i)/N )
%    OffAxis_Resultant_BetweenTrial_SD       = sqrt( sum(n_i*(m_i-gm)^2)/N )
%    OffAxis_Resultant_TotalSD_Check         = sqrt(within^2 + between^2)
%                                               (sanity check -- should
%                                               match OffAxis_Resultant_SD)
%    OffAxis_Resultant_TrialEqualWeighted_Mean = mean(m_i) UNWEIGHTED by
%                                               trial sample count, so a
%                                               longer trial can't
%                                               dominate the mean just by
%                                               having more (correlated)
%                                               cycles -- trial, not
%                                               cycle, is treated as the
%                                               unit of replication here.
%    N_Trials_Pooled_Resamp                  = how many trials went in.
%  OffAxis_Resultant_SD/_Mean (the totals) are still reported and remain
%  valid summaries in their own right -- the within/between split doesn't
%  replace them, it explains what they're made of.
%
%  USAGE
%    KUKA_head_kinematics_summary_AllParticipants(1:10)
%    filt = struct('Condition', 'EO');
%    KUKA_head_kinematics_summary_AllParticipants(1:10, '/data/KUKA', filt)
%    opts = struct('amp_bin_tol_dps', 0.25, 'freq_nominal_targets', [0.4 0.55 0.8 1.0 1.3]);
%    KUKA_head_kinematics_summary_AllParticipants(1:10, [], [], opts)
% ========================================================================

%% ------------------------- CONFIG / DEFAULTS -----------------------------
if nargin < 1 || isempty(participant_nums)
    participant_nums = 1:10;  % match KUKA_analysis_preprocess_driver.m CONFIG range
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
    out_dir = fullfile(base_dir, 'AllParticipants_HeadKinematics');
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ------------------------- PER-PARTICIPANT: DEDUP + POOL + COMPUTE -------
diagnostic_rows = {};
pooled_rows = {};
participants_found  = [];
participants_missing = [];

for pi = 1:numel(participant_nums)
    p = participant_nums(pi);
    csv_out_dir = participant_folder(base_dir, p);

    if ~exist(csv_out_dir, 'dir')
        fprintf('[SKIP participant %d] Folder not found: %s\n', p, csv_out_dir);
        participants_missing(end+1) = p; %#ok<AGROW>
        continue;
    end
    mat_files = dir(fullfile(csv_out_dir, 'trimmed_*.mat'));
    if isempty(mat_files)
        fprintf('[SKIP participant %d] No trimmed_*.mat files in: %s\n', p, csv_out_dir);
        participants_missing(end+1) = p; %#ok<AGROW>
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
        participants_missing(end+1) = p; %#ok<AGROW>
        continue;
    end
    dedup_list = values(best_by_trial);

    fprintf('\n=== Participant %d: %d physical trial(s) ===\n', p, numel(dedup_list));

    % --- Load every deduplicated trial once; build the diagnostic row AND
    % stash the loaded data by (Direction,ParameterSet,Condition) group. ---
    group_key_of = @(tok) strjoin({tok.direction, tok.paramset, tok.condition}, '|');
    group_data  = containers.Map('KeyType', 'char', 'ValueType', 'any');
    group_files = containers.Map('KeyType', 'char', 'ValueType', 'any');
    group_tok   = containers.Map('KeyType', 'char', 'ValueType', 'any');

    n_diag_before = numel(diagnostic_rows);

    for k = 1:numel(dedup_list)
        fpath = fullfile(dedup_list{k}.folder, dedup_list{k}.name);
        tok = dedup_list{k}.tok;

        try
            S = load(fpath, 'data');
            data = S.data;
        catch ME
            fprintf('  [SKIP] Could not load %s: %s\n', dedup_list{k}.name, ME.message);
            continue;
        end
        if ~isfield(data, 'status') || ~strcmp(data.status, 'ok')
            st = 'unknown';
            if isfield(data, 'status'); st = data.status; end
            fprintf('  [SKIP] %s has status "%s" (not ok)\n', dedup_list{k}.name, st);
            continue;
        end
        if ~isfield(data, 'imu_flange') || ~isfield(data, 'imu_head')
            fprintf('  [SKIP] %s is missing imu_flange/imu_head fields.\n', dedup_list{k}.name);
            continue;
        end

        try
            diag_row = compute_trial_kinematics(data, tok, fpath);
            diag_row.Participant = p;
            diagnostic_rows{end+1} = diag_row; %#ok<AGROW>
        catch ME
            fprintf('  [SKIP diagnostic] %s failed: %s\n', dedup_list{k}.name, ME.message);
        end

        gk = group_key_of(tok);
        if ~isKey(group_data, gk)
            group_data(gk) = {};
            group_files(gk) = {};
            group_tok(gk) = tok;
        end
        tmp = group_data(gk); tmp{end+1} = data; group_data(gk) = tmp; %#ok<AGROW>
        tmpf = group_files(gk); tmpf{end+1} = tok.imu_base; group_files(gk) = tmpf; %#ok<AGROW>

        if opts.make_axis_vs_flange_plots
            [t_resamp, ~] = get_resampled_time_base(data);
            plot_trial_axis_vs_flange(t_resamp, ...
                data.imu_head.six_axis_resamp, data.imu_flange.six_axis_resamp, ...
                tok, p, out_dir, 'Resamp');
        end
    end

    n_diag_added = numel(diagnostic_rows) - n_diag_before;
    if n_diag_added > 0
        participants_found(end+1) = p; %#ok<AGROW>
    else
        participants_missing(end+1) = p; %#ok<AGROW>
    end

    % --- Now pool each group and compute the metrics ONCE per group. ---
    gkeys = keys(group_data);
    for gi = 1:numel(gkeys)
        gk = gkeys{gi};
        data_list = group_data(gk);
        files_included = group_files(gk);
        tok0 = group_tok(gk);

        try
            prow = pool_and_compute_group(data_list, tok0);
        catch ME
            fprintf('  [SKIP group] %s failed pooled computation: %s\n', gk, ME.message);
            continue;
        end
        prow.Participant = p;
        prow.N_Trials_Pooled = numel(data_list);
        prow.TrialFiles_Included = strjoin(files_included, ',');
        pooled_rows{end+1} = prow; %#ok<AGROW>

        fprintf('  %-30s : %d trial(s) pooled (%d Resamp samples)\n', gk, numel(data_list), prow.N_Samples_Pooled_Resamp);
    end
end

if isempty(pooled_rows)
    error(['No usable (Participant, Direction, ParameterSet, Condition) groups found under:\n  %s\n' ...
        'Check base_dir, participant_nums, and filt.'], base_dir);
end

%% ------------------------- LEVEL 1: PER-TRIAL DIAGNOSTIC TABLE -----------
% NOT used by anything below -- purely so an individual trial can be
% eyeballed (e.g. an outlier) before trusting the pooled Level A result.
if ~isempty(diagnostic_rows)
    diagnostic_table = struct2table([diagnostic_rows{:}]);
    out_csv_diag = fullfile(out_dir, 'Head_Kinematics_PerTrial_Diagnostic.csv');
    writetable(diagnostic_table, out_csv_diag);
    fprintf('\nLevel 1 (per-trial DIAGNOSTIC ONLY, not used downstream) written to:\n  %s\n', out_csv_diag);
end

pooled_table = struct2table([pooled_rows{:}]);
fprintf('\nComputed pooled kinematics for %d (Participant, Direction, ParameterSet, Condition) group(s) across %d participant(s) (missing/empty: %s).\n', ...
    height(pooled_table), numel(participants_found), mat2str(participants_missing));

%% ------------------------- LEVEL A: PER-PARTICIPANT POOLED TABLE ---------
out_csv_A = fullfile(out_dir, 'Head_Kinematics_PerParticipant_Pooled.csv');
writetable(pooled_table, out_csv_A);
fprintf('Level A (per-participant, pooled-trial kinematics) written to:\n  %s\n', out_csv_A);

%% ------------------- FREQ/AMP BINNING (computed once, used by B, C, D) ----
% AMPLITUDE: tolerance-based clustering (divisive, gap-based) --
% opts.amp_bin_tol_dps, validated against the real data's natural gaps
% (confirmed 0.2 cleanly separates the true amplitude levels).
%
% FREQUENCY: snapped to the known NOMINAL DESIGNED frequencies
% (opts.freq_nominal_targets) instead of inferred from gaps -- the
% frequency values are messier (no clean universal gap size the way
% amplitude has), but the actual designed targets are known, so every
% observed Freq_Hz is assigned to whichever nominal target it's
% numerically closest to. This sidesteps the tolerance-clustering
% ambiguity entirely (no risk of a close pair like 0.38/0.41 landing in
% different bins by accident of which internal gap happened to be
% largest) -- assignment is a straight nearest-target lookup.
%
% Computed on pooled_table BEFORE Level B so that Level B's exact-combo
% table ALSO carries the binned AmpBin/FreqBin columns (used by the
% Level B "MagnitudeByFreqAmp"/"VariabilityByFreqAmp" plots below) --
% previously those plots grouped lines by the raw, unbinned Amp_dps, so
% every near-duplicate designed amplitude (3.93, 4.00, 4.02, 4.08, ...)
% got its own separate legend entry/line instead of being consolidated.
amp_tol = opts.amp_bin_tol_dps;
[amp_bin_id, amp_bin_centers] = cluster_values_tolerance(pooled_table.Amp_dps, amp_tol);
pooled_table.AmpBin = round(amp_bin_centers(amp_bin_id), 3);

if isempty(opts.freq_nominal_targets)
    freq_tol = opts.freq_bin_tol_hz;
    [freq_bin_id, freq_bin_centers] = cluster_values_tolerance(pooled_table.Freq_Hz, freq_tol);
    pooled_table.FreqBin = round(freq_bin_centers(freq_bin_id), 4);
else
    nominal_targets = opts.freq_nominal_targets(:);   % ensure column vector
    freq_bin_id = snap_to_nominal(pooled_table.Freq_Hz, nominal_targets);
    pooled_table.FreqBin = nominal_targets(freq_bin_id);
end

%% ------------------------- LEVEL B: GROUP SUMMARY (exact freq/amp table) --
lb_keys = strcat(pooled_table.ParameterSet, "_", pooled_table.Condition, "_", pooled_table.Direction);
[lb_unique, ~, lb_idx] = unique(lb_keys);

is_numeric_var2 = varfun(@isnumeric, pooled_table, 'OutputFormat', 'uniform');
numeric_vars2   = pooled_table.Properties.VariableNames(is_numeric_var2);
numeric_vars2   = setdiff(numeric_vars2, {'Participant', 'FreqBin', 'AmpBin'}, 'stable');

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
    srow.FreqBin        = sub.FreqBin(1);   % binned label -- used for plot grouping
    srow.AmpBin         = sub.AmpBin(1);    % binned label -- used for plot grouping
    srow.TaskAxes       = sub.TaskAxes(1);
    srow.OffAxisAxes    = sub.OffAxisAxes(1);
    srow.N_Participants = height(sub);
    srow.Participants_Included = strjoin(string(unique(sub.Participant)), ',');

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
    lb_rows{end+1} = srow; %#ok<AGROW>
end
lb_table = struct2table([lb_rows{:}]);
lb_table = sortrows(lb_table, {'Direction', 'Freq_Hz', 'Amp_dps', 'Condition'});

out_csv_B = fullfile(out_dir, 'Head_Kinematics_Group_Summary.csv');
writetable(lb_table, out_csv_B);
fprintf('Level B (group freq/amp summary across participants, exact combo) written to:\n  %s\n', out_csv_B);
disp(lb_table(:, {'Direction','Freq_Hz','Amp_dps','Condition','N_Participants'}));

%% ------------------------- LEVEL C: FREQ x AMP GRID -----------------------
% AmpBin/FreqBin were already computed above (before Level B).
grid_keys = strcat(pooled_table.Direction, "_", pooled_table.Condition, "_", ...
    string(pooled_table.FreqBin), "_", string(pooled_table.AmpBin));
[grid_unique, ~, grid_idx] = unique(grid_keys);

is_numeric_var3 = varfun(@isnumeric, pooled_table, 'OutputFormat', 'uniform');
numeric_vars3   = pooled_table.Properties.VariableNames(is_numeric_var3);
numeric_vars3   = setdiff(numeric_vars3, {'Participant', 'FreqBin', 'AmpBin'}, 'stable');

grid_rows = {};
for g = 1:numel(grid_unique)
    idx = (grid_idx == g);
    sub = pooled_table(idx, :);

    srow = struct();
    srow.Direction   = sub.Direction(1);
    srow.Condition    = sub.Condition(1);
    srow.FreqBin       = sub.FreqBin(1);
    srow.AmpBin          = sub.AmpBin(1);
    srow.TaskAxes          = sub.TaskAxes(1);
    srow.OffAxisAxes         = sub.OffAxisAxes(1);
    srow.N_Participants        = height(sub);
    srow.Participants_Included  = strjoin(string(unique(sub.Participant)), ',');
    srow.ParameterSets_Included  = strjoin(unique(sub.ParameterSet), ',');

    for v = 1:numel(numeric_vars3)
        vn = numeric_vars3{v};
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

out_csv_freq = fullfile(out_dir, 'Head_Kinematics_FreqEffect_at_ConstAmplitude.csv');
writetable(sortrows(grid_table, {'Direction', 'Condition', 'AmpBin', 'FreqBin'}), out_csv_freq);
fprintf('Frequency-effect-at-constant-amplitude table written to:\n  %s\n', out_csv_freq);

out_csv_amp = fullfile(out_dir, 'Head_Kinematics_AmpEffect_at_ConstFrequency.csv');
writetable(sortrows(grid_table, {'Direction', 'Condition', 'FreqBin', 'AmpBin'}), out_csv_amp);
fprintf('Amplitude-effect-at-constant-frequency table written to:\n  %s\n', out_csv_amp);

%% ------------------- SHARED-ACROSS-DIRECTIONS FLAG -----------------------
% Different directions couldn't always achieve the same designed
% frequency/amplitude, so some directions substituted a compromise value
% -- meaning a (FreqBin, AmpBin) cell that exists for one direction often
% doesn't exist at all for another (e.g. AP tested ~3.6/~6.9 dps that no
% other direction ever tested; the diagonals tested ~0.48 Hz that AP/ML
% never tested). A cell is marked SharedAcrossDirections=true only if it
% has data for EVERY direction present in that Condition's data -- this
% is what any fair cross-direction comparison should filter on, not just
% "does this FreqBin/AmpBin label exist somewhere."
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

%% ------------------------- LEVEL D: DIRECTION/CONDITION OVERALL SUMMARY ---
dc_keys = strcat(grid_table.Direction, "_", grid_table.Condition);
[dc_unique, ~, dc_idx] = unique(dc_keys);

metric_cols = { ...
    'OffAxis_to_Task_RMS_Ratio_Resamp_mean', 'Off-axis/task RMS ratio (RECOMMENDED magnitude measure, normalized)'; ...
    'OffAxis_Resultant_RMS_Resamp_mean',     'Off-axis resultant RMS, raw (deg/s)'; ...
    'OffAxis_Resultant_CV_Resamp_mean',      'Off-axis resultant CV (within-trial variability, mean across cells)'};

d4_rows = {};
for g = 1:numel(dc_unique)
    idx = (dc_idx == g);
    sub = grid_table(idx, :);

    drow = struct();
    drow.Direction       = sub.Direction(1);
    drow.Condition        = sub.Condition(1);
    drow.N_FreqAmpCells    = height(sub);

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

out_csv_D = fullfile(out_dir, 'Head_Kinematics_DirectionOverall_Summary.csv');
writetable(direction_summary_table, out_csv_D);
fprintf(['Level D (direction/condition overall summary, irrespective of input) written to:\n  %s\n' ...
    '  NOTE: each direction here is averaged over ITS OWN tested Freq x Amp space, which\n' ...
    '  differs by direction (some directions substituted a compromise value where the\n' ...
    '  designed input could not be achieved) -- so these numbers are NOT directly\n' ...
    '  comparable direction-to-direction. Use Head_Kinematics_DirectionComparison_\n' ...
    '  SharedBinsOnly.csv (below) for that.\n'], out_csv_D);
disp(direction_summary_table);

%% --------- LEVEL D-SHARED: DIRECTION COMPARISON (SHARED BINS ONLY) --------
% Restricts the same grand-mean/AcrossInput-SD/CV calculation to only the
% (FreqBin, AmpBin) cells flagged SharedAcrossDirections=true above --
% i.e. cells every tested direction actually has data for -- so THESE
% numbers ARE a fair like-for-like comparison across directions.
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
        drow.Direction            = sub.Direction(1);
        drow.Condition             = sub.Condition(1);
        drow.N_SharedFreqAmpCells   = height(sub);
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

    out_csv_Dshared = fullfile(out_dir, 'Head_Kinematics_DirectionComparison_SharedBinsOnly.csv');
    writetable(direction_comparison_table, out_csv_Dshared);
    fprintf('Level D-Shared (fair cross-direction comparison, shared Freq/Amp cells only) written to:\n  %s\n', out_csv_Dshared);
    disp(direction_comparison_table);
end

%% ------------- DIAGNOSTIC: (FreqBin, AmpBin) CELLS NOT SHARED --------------
excluded_grid = grid_table(~grid_table.SharedAcrossDirections, :);
if ~isempty(excluded_grid)
    out_csv_excluded = fullfile(out_dir, 'Head_Kinematics_NonSharedFreqAmpCells.csv');
    writetable(sortrows(excluded_grid(:, {'Direction', 'Condition', 'FreqBin', 'AmpBin', 'N_Participants'}), ...
        {'Condition', 'FreqBin', 'AmpBin', 'Direction'}), out_csv_excluded);
    fprintf(['\n%d (FreqBin, AmpBin) cell(s) are NOT shared across every tested direction\n' ...
        '(excluded from Level D-Shared) -- see:\n  %s\n'], height(excluded_grid), out_csv_excluded);
end

%% ------------------------- PLOTS ------------------------------------------
plot_offaxis_freq_amp_effects(lb_table, out_dir);
plot_offaxis_effect(grid_table, out_dir, 'FreqBin', 'AmpBin', ...
    'Input frequency (Hz)', 'Amp %.2f dps', 'OffAxis_FreqEffect_at_ConstAmplitude');
plot_offaxis_effect(grid_table, out_dir, 'AmpBin', 'FreqBin', ...
    'Input amplitude (deg/s)', 'Freq %.3f Hz', 'OffAxis_AmpEffect_at_ConstFrequency');

end % KUKA_head_kinematics_summary_AllParticipants


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function csv_out_dir = participant_folder(base_dir, p)
% Matches KUKA_analysis_preprocess_driver.m's
%   outfile = sprintf('P%02d_sync_summary_order', participant_num);
%   csv_out_dir = fullfile(pwd, outfile);
% If your layout differs, edit this one function.
    csv_out_dir = fullfile(base_dir, sprintf('P%02d_sync_summary_order', p));
end

function filt = fill_default_filter(filt)
    if ~isfield(filt, 'ParameterSet'); filt.ParameterSet = ''; end
    if ~isfield(filt, 'Condition');    filt.Condition    = ''; end
    if ~isfield(filt, 'Direction');    filt.Direction    = ''; end
    if ~isfield(filt, 'Muscle');       filt.Muscle       = []; end
end

function opts = fill_default_opts(opts)
    if ~isfield(opts, 'make_axis_vs_flange_plots')
        opts.make_axis_vs_flange_plots = false;  % batch default: off (see v4 for why)
    end
    if ~isfield(opts, 'save_dir')
        opts.save_dir = '';
    end
    if ~isfield(opts, 'freq_bin_tol_hz')
        opts.freq_bin_tol_hz = 0.02;   % kept for backward compatibility; no
                                        % longer used for FreqBin (see
                                        % freq_nominal_targets below) unless
                                        % you clear freq_nominal_targets
    end
    if ~isfield(opts, 'amp_bin_tol_dps')
        opts.amp_bin_tol_dps = 0.2;    % Level C: +/- deg/s clustered into one AmpBin
    end
    if ~isfield(opts, 'freq_nominal_targets')
        % Known DESIGNED frequency levels (Hz) -- every observed Freq_Hz is
        % snapped to whichever of these it's numerically closest to, rather
        % than inferred from gaps in the data (frequency's real jitter
        % doesn't have one clean universal gap size the way amplitude
        % does). Edit this list if your design changes; set to [] to fall
        % back to tolerance-based clustering via freq_bin_tol_hz instead.
        opts.freq_nominal_targets = [0.4, 0.48, 0.55, 0.8, 1.0, 1.3];
        % 0.48 kept separate from 0.55: 0.48 turned out to be a
        % diagonal-only (DiaL/DiaR) compromise value, distinct from
        % ML's true 0.55 Hz -- see the SharedAcrossDirections logic
        % below, which is what actually matters for cross-direction
        % comparisons, not just how finely frequency gets binned.
    end
end

function bin_id = snap_to_nominal(vals, nominal_targets)
% Assigns each value in vals to whichever entry in nominal_targets it is
% numerically closest to (ties broken toward the lower target). Returns
% an index into nominal_targets per value -- nominal_targets(bin_id) is
% the assigned bin's label. No tolerance/gap logic involved: this is a
% direct nearest-neighbor lookup against KNOWN designed values, used
% when the true experimental setpoints are known (see
% opts.freq_nominal_targets) rather than inferred from the data itself.
    vals = vals(:);
    nominal_targets = nominal_targets(:);
    diffs = abs(vals - nominal_targets');   % N x K matrix of |val - target|
    [~, bin_id] = min(diffs, [], 2);
end

function [bin_id, bin_center] = cluster_values_tolerance(vals, tol)
% DIVISIVE gap-based clustering: looks at the WHOLE sorted set of values
% at once, rather than walking through them once with a single fixed
% anchor or a chain of pairwise comparisons. Starts with everything in
% one bin; while any bin's width (max - min) exceeds tol, it is split at
% its OWN LARGEST INTERNAL GAP (not at a fixed offset from wherever a
% bin happened to start) -- repeated until every bin's width is <= tol.
%
% This directly fixes two different failure modes seen in earlier
% versions of this function:
%  - Anchored-at-bin-start version: a bin's cutoff was bin_start + tol,
%    so two values within tol of EACH OTHER (e.g. 4.02 and 4.09 at
%    tol=0.2) could still be split apart if the bin happened to be
%    anchored at some lower, unrelated value (e.g. 3.85, cutoff 4.05)
%    that put 4.02 in and pushed 4.09 into a new bin.
%  - Chained-to-previous-value version: guaranteed any two values within
%    tol of each other merge, but let a long dense run of close values
%    chain into one bin WIDER than tol.
% Splitting at the largest internal gap avoids both: every bin is
% guaranteed <= tol wide (no unbounded chaining), and a close pair like
% 4.02/4.09 only gets separated if that's genuinely the best (or only)
% place within their group to make a valid split -- not an artifact of
% where the scan happened to start.
    vals = vals(:);
    n = numel(vals);
    [sorted_vals, sort_idx] = sort(vals);

    tol_eps = 1e-9;   % guards against float rounding spuriously splitting a
                      % group whose true span sits exactly at tol (e.g.
                      % 0.39-0.37 not landing on a clean 0.02 in float math)
    clusters = {[1, n]};   % each cluster is a [lo, hi] index range into sorted_vals
    changed = true;
    while changed
        changed = false;
        next_clusters = {};
        for c = 1:numel(clusters)
            lo = clusters{c}(1);
            hi = clusters{c}(2);
            if hi > lo && (sorted_vals(hi) - sorted_vals(lo)) > tol + tol_eps
                seg_gaps = diff(sorted_vals(lo:hi));
                [~, gi] = max(seg_gaps);          % largest internal gap
                split_at = lo + gi - 1;            % last index of the left half
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

function tok = parse_trimmed_filename(fname)
% Parses trimmed_<imu_base>_Muscle<N>_<direction>_<ParameterSet>_<condition>.mat
% back into its component input-parameter fields. Returns [] if fname
% doesn't match.
%  ParameterSet = "<freq*1000>_<amp*100>", e.g. "1800_0403" ->
%    freq_hz = 1.800 Hz, amp_raw = 403, amp_dps = 4.03 deg/s
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
    if ~isempty(filt.Muscle) && tok.muscle ~= filt.Muscle
        tf = false; return;
    end
end

function ref = get_flange_reference(data, direction, six_axis_field)
    rot_field = 'rotation_native';
    if strcmp(six_axis_field, 'six_axis_resamp')
        rot_field = 'rotation_resamp';
    end
    if isfield(data.imu_flange, rot_field) && ~isempty(data.imu_flange.(rot_field))
        ref = data.imu_flange.(rot_field)(:);
    else
        ref = reconstruct_flange_reference(data.imu_flange.(six_axis_field), direction);
    end
end

function ref = reconstruct_flange_reference(flange6, direction)
    d = upper(direction);
    gx = flange6(:, 4);
    gy = flange6(:, 5);

    if strcmp(d, 'ML')
        ref = gx;
    elseif strcmp(d, 'AP')
        ref = gy;
    elseif contains(d, 'DIA')
        gxy = [gx, gy];
        C = cov(gxy);
        [V, D] = eig(C);
        [~, idx] = max(diag(D));
        pc1 = V(:, idx);
        if contains(d, 'R')
            ref_dir = [1; 1];
        else
            ref_dir = [1; -1];
        end
        ref_dir = ref_dir / norm(ref_dir);
        if dot(pc1, ref_dir) < 0
            pc1 = -pc1;
        end
        ref = gxy * pc1;
    else
        ref = gx;
    end
end

function [task_idx, offaxis_idx] = axis_roles_for_direction(direction)
    d = upper(direction);
    if strcmp(d, 'ML')
        task_idx    = 1;
        offaxis_idx = [2 3];
    elseif strcmp(d, 'AP')
        task_idx    = 2;
        offaxis_idx = [1 3];
    elseif contains(d, 'DIA')
        task_idx    = [1 2];
        offaxis_idx = 3;
    else
        task_idx    = [1 2 3];
        offaxis_idx = [];
    end
end

function [peak_corr, peak_lag_sec] = axis_phase_relation(sig, ref, fs)
    peak_corr = NaN;
    peak_lag_sec = NaN;
    sig = sig(:); ref = ref(:);
    n = min(numel(sig), numel(ref));
    if n < 4
        return;
    end
    sig = sig(1:n) - mean(sig(1:n), 'omitnan');
    ref = ref(1:n) - mean(ref(1:n), 'omitnan');
    if all(sig == 0) || all(ref == 0) || any(isnan(sig)) || any(isnan(ref))
        return;
    end
    max_lag = max(1, min(round(n/2), round(2 * fs)));
    try
        [c, lags] = xcorr(sig, ref, max_lag, 'coeff');
    catch
        return;
    end
    [~, idx] = max(abs(c));
    peak_corr    = c(idx);
    peak_lag_sec = lags(idx) / fs;
end

function [t_resamp, fs_resamp] = get_resampled_time_base(data)
    N = size(data.imu_head.six_axis_resamp, 1);
    if isfield(data, 'emg') && isfield(data.emg, 'fs') && ~isempty(data.emg.fs) && isfinite(data.emg.fs)
        fs_resamp = data.emg.fs;
    else
        fs_resamp = 2000;
    end
    t_resamp = (1:N)' / fs_resamp;
end

function m = compute_single_axis_metrics(sig, t, ref, fs)
    sig = sig(:);
    m = struct();
    m.RMS  = sqrt(mean(sig.^2, 'omitnan'));
    m.PtP  = max(sig) - min(sig);
    m.Bias = mean(sig, 'omitnan');
    if numel(t) == numel(sig) && numel(t) > 1
        m.TotalExcursion = trapz(t, abs(sig));
    else
        m.TotalExcursion = NaN;
    end
    [m.PhaseCorr, m.PhaseLag_sec] = axis_phase_relation(sig, ref, fs);
end

function d = compute_offaxis_deviation_metrics(head_gyro, offaxis_idx, t)
    d = struct();
    if isempty(offaxis_idx)
        d.OffAxis_Resultant_RMS            = NaN;
        d.OffAxis_Resultant_Mean           = NaN;
        d.OffAxis_Resultant_PtP            = NaN;
        d.OffAxis_Resultant_TotalExcursion = NaN;
        d.OffAxis_Resultant_SD             = NaN;
        d.OffAxis_Resultant_CV             = NaN;
        return;
    end

    resultant = sqrt(sum(head_gyro(:, offaxis_idx).^2, 2));

    d.OffAxis_Resultant_RMS  = sqrt(mean(resultant.^2, 'omitnan'));
    d.OffAxis_Resultant_Mean = mean(resultant, 'omitnan');
    d.OffAxis_Resultant_PtP  = max(resultant) - min(resultant);
    if numel(t) == numel(resultant) && numel(t) > 1
        d.OffAxis_Resultant_TotalExcursion = trapz(t, resultant);
    else
        d.OffAxis_Resultant_TotalExcursion = NaN;
    end
    d.OffAxis_Resultant_SD = std(resultant, 'omitnan');
    if d.OffAxis_Resultant_Mean ~= 0
        d.OffAxis_Resultant_CV = d.OffAxis_Resultant_SD / d.OffAxis_Resultant_Mean;
    else
        d.OffAxis_Resultant_CV = NaN;
    end
end

function decomp = compute_offaxis_within_between(data_list, offaxis_idx, six_axis_field)
% Within-trial vs between-trial decomposition of the OffAxisResultant
% signal across the trials in data_list (a cell array of per-trial
% `data` structs already filtered to whichever ones contribute to this
% block -- e.g. all of them for Resamp, only the HP-capable ones for HP).
% See file header "WITHIN-TRIAL vs BETWEEN-TRIAL DECOMPOSITION" for the
% formulas. Returns NaNs if there's nothing usable to decompose.
    K = numel(data_list);
    trial_means = nan(K, 1);
    trial_vars  = nan(K, 1);
    trial_n     = nan(K, 1);

    if ~isempty(offaxis_idx)
        for i = 1:K
            d = data_list{i};
            if ~isfield(d, 'imu_head') || ~isfield(d.imu_head, six_axis_field)
                continue;
            end
            gyro = d.imu_head.(six_axis_field)(:, 4:6);
            resultant = sqrt(sum(gyro(:, offaxis_idx).^2, 2));
            resultant = resultant(~isnan(resultant));
            if isempty(resultant)
                continue;
            end
            trial_means(i) = mean(resultant);
            trial_vars(i)  = var(resultant);   % default (N-1) normalization, matches std() used elsewhere
            trial_n(i)     = numel(resultant);
        end
    end

    valid = ~isnan(trial_means) & ~isnan(trial_vars) & trial_n > 0;
    decomp = struct();
    decomp.N_Trials = sum(valid);

    if decomp.N_Trials == 0
        decomp.TrialEqualWeighted_Mean = NaN;
        decomp.WithinTrial_SD          = NaN;
        decomp.BetweenTrial_SD         = NaN;
        decomp.TotalSD_Check           = NaN;
        return;
    end

    tm = trial_means(valid);
    tv = trial_vars(valid);
    tn = trial_n(valid);
    N  = sum(tn);

    grand_mean_weighted = sum(tn .* tm) / N;   % == the concatenated/pooled mean
    within_var  = sum(tn .* tv) / N;

    decomp.TrialEqualWeighted_Mean = mean(tm);           % unweighted across trials
    decomp.WithinTrial_SD          = sqrt(within_var);   % still valid with 1 trial -- it's just that trial's own SD

    if decomp.N_Trials < 2
        % With only one trial, there is nothing to compare it AGAINST --
        % between-trial variability is not zero, it's UNMEASURED. Reporting
        % 0 here would misleadingly claim "perfectly consistent across
        % trials" when the truth is "no trial-to-trial information at
        % all." NaN correctly flags this cell as not estimable, so it
        % doesn't get silently treated as evidence of good consistency.
        decomp.BetweenTrial_SD = NaN;
        decomp.TotalSD_Check   = NaN;
    else
        between_var = sum(tn .* (tm - grand_mean_weighted).^2) / N;
        decomp.BetweenTrial_SD = sqrt(between_var);
        decomp.TotalSD_Check   = sqrt(within_var + between_var);
    end
end

function block = compute_kinematics_block(head6, flange6, ref, fs, t, task_idx, offaxis_idx, suffix)
    accel_names = {'Ax', 'Ay', 'Az'};
    axis_names  = {'Gx', 'Gy', 'Gz'};
    block = struct();

    for c = 1:3
        block.(['Flange_' accel_names{c} '_mean' suffix]) = mean(flange6(:, c), 'omitnan');
        block.(['Flange_' accel_names{c} '_std' suffix])  = std(flange6(:, c), 'omitnan');
    end
    for c = 1:3
        block.(['Flange_' axis_names{c} '_mean' suffix]) = mean(flange6(:, 3 + c), 'omitnan');
        block.(['Flange_' axis_names{c} '_std' suffix])  = std(flange6(:, 3 + c), 'omitnan');
    end
    for c = 1:3
        block.(['Head_' accel_names{c} '_mean' suffix]) = mean(head6(:, c), 'omitnan');
        block.(['Head_' accel_names{c} '_std' suffix])  = std(head6(:, c), 'omitnan');
    end

    head_gyro = head6(:, 4:6);
    axis_rms = nan(1, 3);
    for c = 1:3
        m = compute_single_axis_metrics(head_gyro(:, c), t, ref, fs);
        fn = axis_names{c};
        block.(['Head_' fn '_RMS' suffix])            = m.RMS;
        block.(['Head_' fn '_PtP' suffix])            = m.PtP;
        block.(['Head_' fn '_TotalExcursion' suffix]) = m.TotalExcursion;
        block.(['Head_' fn '_Bias' suffix])           = m.Bias;
        block.(['Head_' fn '_PhaseCorr' suffix])      = m.PhaseCorr;
        block.(['Head_' fn '_PhaseLag_sec' suffix])   = m.PhaseLag_sec;
        axis_rms(c) = m.RMS;
    end

    if ~isempty(task_idx) && ~isempty(offaxis_idx)
        block.(['OffAxis_to_Task_RMS_Ratio' suffix]) = mean(axis_rms(offaxis_idx), 'omitnan') / mean(axis_rms(task_idx), 'omitnan');
    else
        block.(['OffAxis_to_Task_RMS_Ratio' suffix]) = NaN;
    end

    % --- Per-calibrated-axis gyro ENERGY FRACTION (unsigned, magnitude-
    % based -- the RMS-script counterpart to the PCA script's directional
    % DeviationVec_Gx/Gy/Gz). Head_Gx/Gy/Gz_RMS above already give the
    % per-axis magnitude individually, but nothing yet shows what
    % FRACTION of the head's total rotational energy landed on each axis
    % -- this is direction-agnostic (works the same way for ML/AP/
    % diagonal, unlike trying to special-case which columns count as
    % "off-axis" per direction) and always sums to 1 across the three axes.
    total_energy = sum(axis_rms.^2, 'omitnan');
    for c = 1:3
        fn = axis_names{c};
        if total_energy > 0
            block.(['Head_GyroEnergyFrac_' fn suffix]) = (axis_rms(c)^2) / total_energy;
        else
            block.(['Head_GyroEnergyFrac_' fn suffix]) = NaN;
        end
    end

    dev = compute_offaxis_deviation_metrics(head_gyro, offaxis_idx, t);
    dev_fn = fieldnames(dev);
    for i = 1:numel(dev_fn)
        block.([dev_fn{i} suffix]) = dev.(dev_fn{i});
    end

    block.(['Head_Overall_GyroMag_RMS' suffix]) = sqrt(mean(sum(head_gyro.^2, 2), 'omitnan'));
    if numel(t) > 1
        block.(['Duration_sec' suffix]) = t(end) - t(1);
    else
        block.(['Duration_sec' suffix]) = NaN;
    end
end

function row = compute_trial_kinematics(data, tok, fpath)
% SINGLE-TRIAL metrics -- used ONLY to build the Level 1 diagnostic table.
% Every aggregate level uses pool_and_compute_group instead.
    row = struct();
    row.IMU_FileName = string(tok.imu_base);
    row.Muscle       = tok.muscle;
    row.Direction    = string(tok.direction);
    row.ParameterSet = string(tok.paramset);
    row.Condition    = string(tok.condition);
    row.Freq_Hz      = tok.freq_hz;
    row.Amp_raw      = tok.amp_raw;
    row.Amp_dps       = tok.amp_dps;
    row.MatPath      = string(fpath);

    axis_names = {'Gx', 'Gy', 'Gz'};
    [task_idx, offaxis_idx] = axis_roles_for_direction(tok.direction);
    row.TaskAxes    = strjoin(axis_names(task_idx), ',');
    row.OffAxisAxes = strjoin(axis_names(offaxis_idx), ',');

    [t_resamp, fs_resamp] = get_resampled_time_base(data);
    ref_resamp = get_flange_reference(data, tok.direction, 'six_axis_resamp');
    resamp_block = compute_kinematics_block(data.imu_head.six_axis_resamp, ...
        data.imu_flange.six_axis_resamp, ref_resamp, fs_resamp, t_resamp, ...
        task_idx, offaxis_idx, '_Resamp');
    row = merge_struct(row, resamp_block);
end

function row = pool_and_compute_group(data_list, tok0)
% POOLED metrics: concatenates every trial's six_axis_resamp
% head/flange/reference samples end-to-end,
% builds a synthetic monotonic time base over the pooled length (same
% approach each trial's own resampled time base already used -- sample
% count over rate, not a real shared clock, so concatenation across
% trials is no less valid than within one trial), and calls
% compute_kinematics_block ONCE per block on the pooled arrays. This is
% what every aggregate level (A/B/C/D) is actually built from.
    row = struct();
    row.Direction    = string(tok0.direction);
    row.ParameterSet = string(tok0.paramset);
    row.Condition    = string(tok0.condition);
    row.Freq_Hz      = tok0.freq_hz;
    row.Amp_raw      = tok0.amp_raw;
    row.Amp_dps      = tok0.amp_dps;

    axis_names = {'Gx', 'Gy', 'Gz'};
    [task_idx, offaxis_idx] = axis_roles_for_direction(tok0.direction);
    row.TaskAxes    = strjoin(axis_names(task_idx), ',');
    row.OffAxisAxes = strjoin(axis_names(offaxis_idx), ',');

    % --- Resamp block: pool across every trial (all should have this) ---
    head_pool = []; flange_pool = []; ref_pool = [];
    fs_resamp_used = NaN;
    for i = 1:numel(data_list)
        d = data_list{i};
        [~, fs_i] = get_resampled_time_base(d);
        if isnan(fs_resamp_used)
            fs_resamp_used = fs_i;
        end
        ref_i = get_flange_reference(d, tok0.direction, 'six_axis_resamp');
        head_pool   = [head_pool;   d.imu_head.six_axis_resamp];   %#ok<AGROW>
        flange_pool = [flange_pool; d.imu_flange.six_axis_resamp]; %#ok<AGROW>
        ref_pool    = [ref_pool;    ref_i(:)];                     %#ok<AGROW>
    end
    N_resamp = size(head_pool, 1);
    t_resamp_pool = (0:N_resamp-1)' / fs_resamp_used;
    resamp_block = compute_kinematics_block(head_pool, flange_pool, ref_pool, ...
        fs_resamp_used, t_resamp_pool, task_idx, offaxis_idx, '_Resamp');
    row = merge_struct(row, resamp_block);
    row.N_Samples_Pooled_Resamp = N_resamp;

    decomp_resamp = compute_offaxis_within_between(data_list, offaxis_idx, 'six_axis_resamp');
    row.OffAxis_Resultant_TrialEqualWeighted_Mean_Resamp = decomp_resamp.TrialEqualWeighted_Mean;
    row.OffAxis_Resultant_WithinTrial_SD_Resamp           = decomp_resamp.WithinTrial_SD;
    row.OffAxis_Resultant_BetweenTrial_SD_Resamp          = decomp_resamp.BetweenTrial_SD;
    row.OffAxis_Resultant_TotalSD_Check_Resamp            = decomp_resamp.TotalSD_Check;
    row.N_Trials_Pooled_Resamp                             = decomp_resamp.N_Trials;
end

function s = merge_struct(s, extra)
    fn = fieldnames(extra);
    for i = 1:numel(fn)
        s.(fn{i}) = extra.(fn{i});
    end
end

function plot_trial_axis_vs_flange(t, head6, flange6, tok, participant, out_dir, label)
    t = t(:);
    n = min([numel(t), size(head6, 1), size(flange6, 1)]);
    t = t(1:n); head6 = head6(1:n, :); flange6 = flange6(1:n, :);

    [task_idx, ~] = axis_roles_for_direction(tok.direction);
    axis_names = {'Gx', 'Gy', 'Gz'};
    role_label = repmat({'OFF-AXIS'}, 1, 3);
    role_label(task_idx) = {'TASK'};

    plot_dir = fullfile(out_dir, 'AxisVsFlange_Plots', label);
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir);
    end

    fname = sprintf('AxisVsFlange_P%02d_%s_Muscle%d_%s_%s_%s_%s.png', ...
        participant, tok.imu_base, tok.muscle, tok.direction, tok.paramset, tok.condition, label);
    out_png = fullfile(plot_dir, fname);

    try
        fig = figure('Visible', 'off', 'Position', [100 100 1000 800]);
        ax = gobjects(1, 3);
        for c = 1:3
            ax(c) = subplot(3, 1, c);
            plot(t, flange6(:, 3 + c), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'DisplayName', 'Flange');
            hold on;
            plot(t, head6(:, 3 + c), 'Color', [0.10 0.30 0.80], 'LineWidth', 1.0, 'DisplayName', 'Head');
            hold off;
            ylabel(sprintf('%s (deg/s)', axis_names{c}));
            title(sprintf('%s -- %s', axis_names{c}, role_label{c}));
            grid on;
            if c == 1
                legend('Location', 'best');
            end
            if c == 3
                xlabel('Time (s)');
            end
        end
        linkaxes(ax, 'y');
        all_gyro = [head6(:, 4:6); flange6(:, 4:6)];
        max_peak = max(abs(all_gyro(:))) * 1.1;
        ylim(ax(1), [-max_peak, max_peak]);
        sgtitle(sprintf('P%02d | %s | Muscle %d | %s | %s | %s | %s', ...
            participant, tok.imu_base, tok.muscle, tok.direction, tok.paramset, tok.condition, label), ...
            'Interpreter', 'none');
        print(fig, out_png, '-dpng', '-r150');
        close(fig);
    catch ME
        fprintf('[WARNING] Could not generate %s axis-vs-flange plot for P%02d %s: %s\n', ...
            label, participant, tok.imu_base, ME.message);
        if exist('fig', 'var') && isvalid(fig)
            close(fig);
        end
    end
end

function plot_offaxis_freq_amp_effects(lb_table, out_dir)
% One PNG per Condition present in lb_table (EO plotted separately from
% EC, never overlaid on the same axes/line) -- mixing them would make an
% EO point and an EC point at the same frequency look like two noisy
% samples of one line, when they're actually the comparison of interest.
    try
        conditions = unique(lb_table.Condition, 'stable');
        directions = unique(lb_table.Direction, 'stable');
        metrics = { ...
            'OffAxis_Resultant_RMS_Resamp',  'Off-axis resultant RMS (deg/s)',  'OffAxis_MagnitudeByFreqAmp'; ...
            'OffAxis_Resultant_CV_Resamp',   'Off-axis resultant CV (SD/Mean)', 'OffAxis_VariabilityByFreqAmp' };

        for ci = 1:numel(conditions)
            this_cond = conditions(ci);
            ctable = lb_table(lb_table.Condition == this_cond, :);
            if isempty(ctable)
                continue;
            end

            for mi = 1:size(metrics, 1)
                mean_col = [metrics{mi,1} '_mean'];
                sd_col   = [metrics{mi,1} '_between_SD'];
                if ~ismember(mean_col, ctable.Properties.VariableNames)
                    continue;
                end

                fig = figure('Visible', 'off', 'Position', [100 100 1400 350 * numel(directions)]);
                for di = 1:numel(directions)
                    subplot(numel(directions), 1, di);
                    dsub = ctable(ctable.Direction == directions(di), :);
                    amps = unique(dsub.AmpBin);   % binned, not raw Amp_dps -- see header note
                    hold on;
                    for ai = 1:numel(amps)
                        asub = dsub(dsub.AmpBin == amps(ai), :);
                        asub = sortrows(asub, 'Freq_Hz');
                        errorbar(asub.Freq_Hz, asub.(mean_col), asub.(sd_col), '-o', ...
                            'DisplayName', sprintf('Amp %.2f dps', amps(ai)));
                    end
                    hold off;
                    grid on;
                    xlabel('Input frequency (Hz)');
                    ylabel(metrics{mi,2});
                    title(sprintf('Direction: %s', directions(di)));
                    legend('Location', 'best');
                end
                sgtitle(sprintf('%s -- Condition: %s', strrep(metrics{mi,2}, '_', '\_'), this_cond));
                out_png = fullfile(out_dir, sprintf('%s_%s.png', metrics{mi,3}, this_cond));
                print(fig, out_png, '-dpng', '-r150');
                close(fig);
                fprintf('Saved: %s\n', out_png);
            end
        end
    catch ME
        fprintf('[WARNING] Could not generate off-axis freq/amp plots: %s\n', ME.message);
    end
end

function plot_offaxis_effect(grid_table, out_dir, x_field, line_field, x_label, line_label_fmt, file_prefix)
% One PNG per Condition present in grid_table (EO plotted separately
% from EC). Previously this drew every condition's points onto the same
% line for a given AmpBin/FreqBin, since dsub/line grouping never
% checked Condition -- an EO point and an EC point at the same nominal
% amplitude/frequency would land on the same line indistinguishably.
% Fixed to loop over Condition first and write a separate figure/file
% per condition, matching KUKA_head_PCA_deviation_AllParticipants.m's
% plot_pca_deviation_effect.
    try
        conditions = unique(grid_table.Condition, 'stable');
        directions = unique(grid_table.Direction, 'stable');
        metrics = { ...
            'OffAxis_Resultant_RMS_Resamp', 'Off-axis resultant RMS (deg/s)',  [file_prefix '_Magnitude']; ...
            'OffAxis_Resultant_CV_Resamp',  'Off-axis resultant CV (SD/Mean)', [file_prefix '_Variability'] };

        for ci = 1:numel(conditions)
            this_cond = conditions(ci);
            ctable = grid_table(grid_table.Condition == this_cond, :);
            if isempty(ctable)
                continue;
            end

            for mi = 1:size(metrics, 1)
                mean_col = [metrics{mi,1} '_mean'];
                sd_col   = [metrics{mi,1} '_between_SD'];
                if ~ismember(mean_col, ctable.Properties.VariableNames)
                    continue;
                end

                fig = figure('Visible', 'off', 'Position', [100 100 1400 350 * numel(directions)]);
                for di = 1:numel(directions)
                    subplot(numel(directions), 1, di);
                    dsub = ctable(ctable.Direction == directions(di), :);
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
                    title(sprintf('Direction: %s', directions(di)));
                    legend('Location', 'best');
                end
                sgtitle(sprintf('%s -- Condition: %s', strrep(metrics{mi,2}, '_', '\_'), this_cond));
                out_png = fullfile(out_dir, sprintf('%s_%s.png', metrics{mi,3}, this_cond));
                print(fig, out_png, '-dpng', '-r150');
                close(fig);
                fprintf('Saved: %s\n', out_png);
            end
        end
    catch ME
        fprintf('[WARNING] Could not generate %s plots: %s\n', file_prefix, ME.message);
    end
end
