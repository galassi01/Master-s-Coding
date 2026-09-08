function KUKA_head_kinematics_summary_v4(csv_out_dir, filt, opts)
%% ========================================================================
%  HEAD KINEMATICS + OFF-AXIS COMPENSATION ERROR SUMMARY
% ========================================================================
%  Consumes the .mat files written by KUKA_analysis_preprocess_consumer.m
%  (one per participant / IMU-file / Muscle bank, run in a loop by
%  KUKA_analysis_preprocess_driver.m), named:
%
%     trimmed_<imu_base>_Muscle<N>_<direction>_<ParameterSet>_<condition>.mat
%
%  e.g. trimmed_KUKA_neckEMG_P06_Trial_Block1_01_Muscle1_AP_1800_0504_EO.mat
%
%  Each file contains a struct `data` with (among other things):
%     data.imu_flange.six_axis_resamp  [N x 6] = [Ax Ay Az Gx Gy Gz], the
%                                       primary flange data used here --
%                                       data.imu_flange.six_axis_native
%                                       is no longer saved by the consumer
%     data.imu_head.six_axis_native    [N x 6] = [Ax Ay Az Gx Gy Gz]
%                                       (kept in the .mat file, but unused
%                                       by this script -- see below)
%     data.imu_head.t_native, data.imu_head.fs_native
%     data.imu_flange.rotation_native  direction-projected perturbation
%                                       signal delivered to the flange
%                                       (used as the PhaseCorr reference
%                                       when present; reconstructed from
%                                       six_axis_resamp otherwise -- see
%                                       get_flange_reference() below)
%
%  These six-axis arrays are already the "movement-trimmed" data -- i.e.
%  Level 3 in the consumer's terminology: sync-trimmed, then pulled in by
%  movement_buffer_sec (settle period removed), then trimmed to whole
%  physical swing cycles. That is the only level this pipeline version
%  saves to disk, so "movement-trimmed" and "what's in the .mat file" are
%  the same thing here.
%
%  IMPORTANT: ParameterSet and Condition are NOT stored as fields inside
%  data.meta -- the consumer only bakes them into the output filename.
%  This script recovers them (plus Direction and Muscle) by parsing the
%  filename, and carries them through as explicit columns so every result
%  stays linked to the input-parameter condition it came from.
%
%  WHAT THIS SCRIPT QUANTIFIES
%  ----------------------------
%  The task delivers a perturbation around one axis (ML -> roughly Gx,
%  AP -> roughly Gy, diagonal -> a Gx+Gy combination) that requires
%  compensation. The behavioral error of interest is head rotation on the
%  WRONG axis -- i.e. compensating with the wrong movement. Per trial,
%  each head gyro axis (Gx, Gy, Gz) is labeled TASK (intended for this
%  direction) or OFF-AXIS (should not have been recruited), per this
%  mapping:
%
%     Direction     Task axis (axes)      Off-axis (error) axis(es)
%     ML            Gx                    Gy, Gz
%     AP            Gy                    Gx, Gz
%     DiaR / DiaL   Gx  AND  Gy (each      Gz
%                   analyzed separately)
%
%  Gz (yaw) is never a task axis, since ML/AP/diagonal perturbations are
%  all delivered in the horizontal roll/pitch plane.
%
%  For EVERY axis (Gx, Gy, Gz) -- not just the off-axis ones, so the task
%  axis is available as a same-trial baseline -- this script computes:
%    - RMS                 : overall magnitude of rotation on that axis
%    - PtP                 : peak-to-peak range on that axis
%    - TotalExcursion      : cumulative |angular| travel on that axis
%                            (trapz of |signal|; does not cancel like a
%                            plain net/signed integral would)
%    - Bias                : signed mean -- a persistent directional lean
%                            on that axis (e.g. consistently drifting
%                            clockwise rather than swinging symmetrically)
%    - PhaseCorr / PhaseLag: cross-correlation (peak, and the lag in
%                            seconds it occurs at) between that head axis
%                            and the flange's direction-projected
%                            rotation_native signal, i.e. the actual
%                            perturbation delivered to the arm. High
%                            correlation at a small/consistent lag on an
%                            OFF-axis means the error is a coupled,
%                            perturbation-driven compensation mistake
%                            rather than incidental noise.
%
%  Per trial it also reports OffAxis_to_Task_RMS_Ratio: mean off-axis RMS
%  divided by mean task-axis RMS -- a single number for "how much wrong-
%  axis movement relative to correct-axis movement," useful for ranking/
%  comparing across parameter sets.
%
%  UNITS NOTE: accel/gyro units follow whatever apply_mems_calibration /
%  calibrateIMU produce in the consumer (typically g for accel and deg/s
%  for gyro for this sensor family) -- this script does not rescale
%  anything.
%
%  DEPENDENCY: PhaseCorr/PhaseLag use xcorr() (Signal Processing Toolbox).
%  If unavailable, those two columns will come back as NaN (see
%  axis_phase_relation, wrapped in try/catch) but everything else still
%  runs.
%
%  PER-TRIAL AXIS-VS-FLANGE PLOTS: for every trial, TWO 3-panel time-
%  series figures are saved (Gx/Gy/Gz, head overlaid on its matching
%  flange axis, labeled TASK/OFF-AXIS) to AxisVsFlange_Plots/Native and
%  AxisVsFlange_Plots/Resamp under csv_out_dir -- these are the visual
%  counterpart to PhaseCorr/PhaseLag: you can see directly whether a head
%  axis tracks the flange's timing or not. Set opts.make_axis_vs_flange_plots
%  = false to skip these (e.g. for a large batch where you only want the
%  summary CSVs).
%
%  RESAMPLED vs GYRO-HIGH-PASS: every metric above is computed on TWO
%  versions of the six-axis data per trial (a NATIVE, unfiltered-at-
%  native-imu_fs version previously existed as a third, unsuffixed block,
%  but was retired along with data.imu_flange.six_axis_native, which the
%  consumer no longer saves -- data.imu_head.six_axis_native itself is
%  untouched and still in the .mat file, just unused here now):
%    _Resamp      data.imu_*.six_axis_resamp  -- 5 Hz lowpass, resampled
%                 onto the EMG time base by the consumer (Section 8b).
%                 This is now the PRIMARY block for head kinematics.
%    _HP          data.imu_*.six_axis_native_hp -- native imu_fs; accel
%                 unchanged, gyro run through a 0.1 Hz zero-phase
%                 high-pass (consumer Section 3b) to strip the near-DC
%                 residual gyro calibration bias while leaving genuine
%                 task-frequency movement untouched. Applied to the FULL
%                 sync-trimmed signal before movement-buffer/cycle-trim
%                 slicing, so filtfilt's edge transients settle in the
%                 buffer padding rather than landing inside the analysis
%                 window. Only present if the .mat file came from a
%                 consumer version that saves six_axis_native_hp -- older
%                 (or this pipeline's current) files get NaN-filled _HP
%                 columns (row.GyroHP_Available flags which case it was)
%                 rather than erroring, so a mixed-vintage batch still
%                 processes cleanly. Its PhaseCorr/PhaseLag reference is
%                 the flange's OWN six_axis_native_hp (not six_axis_native,
%                 which is retired), so a bias-corrected head signal is
%                 compared against an equivalently bias-corrected flange
%                 reference.
%                 IMPORTANT: Head_G*_Bias_HP is expected to land near
%                 zero BY CONSTRUCTION (that's what a high-pass filter
%                 does) -- it's a sanity check the filter worked, not an
%                 independent finding. The metrics that matter from this
%                 block are RMS/PtP/TotalExcursion_HP.
%
%  Because six_axis_resamp's sample spacing is fixed and not tied to the
%  original recording clock, its time axis is reconstructed as a plain
%  sample-count-over-rate vector, t = (1:N)/fs_resamp -- fs_resamp is
%  taken from data.emg.fs when present (the rate the consumer actually
%  resampled onto) and falls back to 2000 Hz (the driver's default
%  emg_fs) only if that field is missing.
%
%  USAGE
%    % Process every trimmed_*.mat file found (script-mode CONFIG below):
%    KUKA_head_kinematics_summary()
%
%    % Process a specific participant's output folder:
%    KUKA_head_kinematics_summary('/path/to/P10_sync_summary_order')
%
%    % Restrict to one input-parameter condition:
%    filt = struct('ParameterSet', '1800_0504', 'Condition', 'EO');
%    KUKA_head_kinematics_summary('/path/to/P10_sync_summary_order', filt);
%
%    % Restrict further by direction and/or muscle bank if needed:
%    filt = struct('ParameterSet', '1800_0504', 'Condition', 'EO', ...
%                   'Direction', 'AP', 'Muscle', 1);
%
%    % Skip the per-trial axis-vs-flange plots (summary CSVs only):
%    opts = struct('make_axis_vs_flange_plots', false);
%    KUKA_head_kinematics_summary('/path/to/P10_sync_summary_order', [], opts);
% ========================================================================

%% ------------------------- CONFIG (script-mode defaults) ----------------
if nargin < 1 || isempty(csv_out_dir)
    participant_num = 10;  % match KUKA_analysis_preprocess_driver.m CONFIG
    outfile     = sprintf('P%02d_sync_summary_order', participant_num);
    csv_out_dir = fullfile(pwd, outfile);
end

if nargin < 2 || isempty(filt)
    filt = struct();  % empty = no filtering, process every trimmed_*.mat file
end
filt = fill_default_filter(filt);

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = fill_default_opts(opts);

%% ------------------------- FIND TRIMMED DATA FILES -----------------------
mat_files = dir(fullfile(csv_out_dir, 'trimmed_*.mat'));
if isempty(mat_files)
    error(['No trimmed_*.mat files found in:\n  %s\n' ...
        'Run KUKA_analysis_preprocess_driver.m (which calls ' ...
        'KUKA_analysis_preprocess_consumer.m) for this participant first.'], ...
        csv_out_dir);
end
fprintf('Found %d trimmed data file(s) in:\n  %s\n', numel(mat_files), csv_out_dir);

%% ------------------------- PER-TRIAL KINEMATICS --------------------------
rows = {};
for k = 1:numel(mat_files)
    fpath = fullfile(mat_files(k).folder, mat_files(k).name);

    tok = parse_trimmed_filename(mat_files(k).name);
    if isempty(tok)
        fprintf('[SKIP] Could not parse filename (unexpected format): %s\n', mat_files(k).name);
        continue;
    end

    if ~passes_filter(tok, filt)
        continue;
    end

    try
        S = load(fpath, 'data');
        data = S.data;
    catch ME
        fprintf('[SKIP] Could not load %s: %s\n', mat_files(k).name, ME.message);
        continue;
    end

    if ~isfield(data, 'status') || ~strcmp(data.status, 'ok')
        st = 'unknown';
        if isfield(data, 'status'); st = data.status; end
        fprintf('[SKIP] %s has status "%s" (not ok)\n', mat_files(k).name, st);
        continue;
    end

    if ~isfield(data, 'imu_flange') || ~isfield(data, 'imu_head')
        fprintf('[SKIP] %s is missing imu_flange/imu_head fields.\n', mat_files(k).name);
        continue;
    end

    try
        row = compute_trial_kinematics(data, tok, fpath);
        rows{end+1} = row; %#ok<AGROW>
    catch ME
        fprintf('[SKIP] %s failed kinematics computation: %s\n', mat_files(k).name, ME.message);
        continue;
    end

    if opts.make_axis_vs_flange_plots
        % Native plot retired along with data.imu_flange.six_axis_native
        % (no longer saved by the consumer) -- Resamp is now the primary
        % axis-vs-flange visual, GyroHP still available when present.
        [t_resamp, ~] = get_resampled_time_base(data);
        plot_trial_axis_vs_flange(t_resamp, ...
            data.imu_head.six_axis_resamp, data.imu_flange.six_axis_resamp, ...
            tok, csv_out_dir, 'Resamp');

        if isfield(data.imu_head, 'six_axis_native_hp') && isfield(data.imu_flange, 'six_axis_native_hp')
            plot_trial_axis_vs_flange(data.imu_head.t_native(:), ...
                data.imu_head.six_axis_native_hp, data.imu_flange.six_axis_native_hp, ...
                tok, csv_out_dir, 'GyroHP');
        end
    end
end

if isempty(rows)
    error(['No trials matched the requested filter (or none had usable data). ' ...
        'Check csv_out_dir / filt.ParameterSet / filt.Condition / filt.Direction / filt.Muscle.']);
end

trial_table = struct2table([rows{:}]);
fprintf('\nComputed kinematics for %d matching trial(s).\n', height(trial_table));

%% ------------------------- SAVE PER-TRIAL TABLE ---------------------------
out_csv_trials = fullfile(csv_out_dir, 'Head_Kinematics_PerTrial.csv');
writetable(trial_table, out_csv_trials);
fprintf('Per-trial kinematics written to:\n  %s\n', out_csv_trials);

%% ------------------------- AGGREGATE ACROSS TRIALS ------------------------
% Groups by ParameterSet + Condition + Direction. This collapses Muscle
% 1/2 bank duplicates and any repeated trials that share the same input-
% parameter condition into a single "average kinematics/error for this
% condition" row.
group_keys = strcat(trial_table.ParameterSet, "_", trial_table.Condition, "_", trial_table.Direction);
[unique_keys, ~, grp_idx] = unique(group_keys); %#ok<ASGLU>

is_numeric_var = varfun(@isnumeric, trial_table, 'OutputFormat', 'uniform');
numeric_vars   = trial_table.Properties.VariableNames(is_numeric_var);

summary_rows = {};
for g = 1:numel(unique_keys)
    idx = (grp_idx == g);
    sub = trial_table(idx, :);

    srow = struct();
    srow.ParameterSet = sub.ParameterSet(1);
    srow.Condition    = sub.Condition(1);
    srow.Direction    = sub.Direction(1);
    srow.TaskAxes     = sub.TaskAxes(1);
    srow.OffAxisAxes  = sub.OffAxisAxes(1);
    srow.N_Trials     = height(sub);

    for v = 1:numel(numeric_vars)
        vn = numeric_vars{v};
        if strcmp(vn, 'Muscle')
            srow.Muscles_Included = strjoin(string(unique(sub.Muscle)), ',');
            continue;
        end
        srow.(vn) = mean(sub.(vn), 'omitnan');
    end
    summary_rows{end+1} = srow; %#ok<AGROW>
end
summary_table = struct2table([summary_rows{:}]);

out_csv_summary = fullfile(csv_out_dir, 'Head_Kinematics_Summary.csv');
writetable(summary_table, out_csv_summary);
fprintf('Condition-level average kinematics/error written to:\n  %s\n', out_csv_summary);

%% ------------------------- REPORT + PLOTS ---------------------------------
disp(summary_table);
plot_head_rotation_overview(summary_table, csv_out_dir);

end % KUKA_head_kinematics_summary


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function filt = fill_default_filter(filt)
% Fills in any filter fields the caller didn't supply, so downstream code
% can always assume they exist. Empty ('' or []) means "don't filter on
% this field".
    if ~isfield(filt, 'ParameterSet'); filt.ParameterSet = ''; end
    if ~isfield(filt, 'Condition');    filt.Condition    = ''; end
    if ~isfield(filt, 'Direction');    filt.Direction    = ''; end
    if ~isfield(filt, 'Muscle');       filt.Muscle       = []; end
end

function opts = fill_default_opts(opts)
% Fills in any opts fields the caller didn't supply.
    if ~isfield(opts, 'make_axis_vs_flange_plots')
        opts.make_axis_vs_flange_plots = true;
    end
end

function tok = parse_trimmed_filename(fname)
% Parses trimmed_<imu_base>_Muscle<N>_<direction>_<ParameterSet>_<condition>.mat
% (the exact pattern KUKA_analysis_preprocess_consumer.m writes) back into
% its component input-parameter fields. Returns [] if fname doesn't match.
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
    else
        tok.freq_hz = NaN;
        tok.amp_raw = NaN;
    end
end

function tf = passes_filter(tok, filt)
% Case-insensitive exact match against any non-empty filter field.
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
% Returns the flange reference signal used for PhaseCorr/PhaseLag.
% Prefers the consumer's own rotation_native/rotation_resamp field when
% present (best fidelity, matches whichever direction-projection method
% the consumer used, PCA fit on the WIDE sync-trimmed signal for
% diagonal trials). Falls back to reconstructing it from six_axis_native/
% six_axis_resamp when that field has been trimmed from the .mat file --
% see reconstruct_flange_reference() for the ML/AP/diagonal logic and its
% accuracy caveat for diagonal trials.
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
% Rebuilds the direction-projected flange reference signal directly from
% six-axis data, for use when rotation_native/rotation_resamp aren't
% saved in the .mat file. Mirrors the consumer's own projection logic:
%   ML  -> raw Gx
%   AP  -> raw Gy
%   Dia -> PC1 of [Gx, Gy], sign-fixed to match the R/L diagonal label
%
% CAVEAT (diagonal trials only): the consumer fits this PCA on the WIDE
% sync-trimmed signal (many swing cycles). Here it can only be fit on
% whatever six-axis data is available downstream -- for six_axis_native/
% six_axis_resamp that's just the narrow, already cycle-trimmed movement
% segment (often close to a single cycle). Less data going into the PCA
% fit means a noisier estimate of the true swing-plane orientation than
% the consumer's own rotation_native/rotation_resamp -- a real accuracy
% gap, not just a formula difference. ML/AP reconstruction doesn't have
% this problem (no PCA involved, just a raw column).
%
% PCA is done manually via eig(cov(...)) rather than MATLAB's pca() to
% avoid a Statistics and Machine Learning Toolbox dependency -- for 2
% columns this gives the identical PC1 (up to sign, which is fixed below
% exactly as the consumer does).
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
            ref_dir = [1; 1];    % forward + right
        else
            ref_dir = [1; -1];   % forward + left
        end
        ref_dir = ref_dir / norm(ref_dir);
        if dot(pc1, ref_dir) < 0
            pc1 = -pc1;
        end
        ref = gxy * pc1;
    else
        % Unrecognized direction label -- fall back to Gx rather than guess.
        ref = gx;
    end
end

function [task_idx, offaxis_idx] = axis_roles_for_direction(direction)
% Maps a trial's perturbation direction onto which head gyro axis/axes
% (1=Gx, 2=Gy, 3=Gz) are TASK (intended) vs OFF-AXIS (compensation
% error). Gz (yaw) is never a task axis, since ML/AP/diagonal
% perturbations are all delivered in the horizontal roll/pitch plane.
% Diagonal trials treat Gx and Gy as two SEPARATE task axes (not a single
% PCA-collapsed axis), per design decision.
    d = upper(direction);
    if strcmp(d, 'ML')
        task_idx    = 1;       % Gx
        offaxis_idx = [2 3];   % Gy, Gz
    elseif strcmp(d, 'AP')
        task_idx    = 2;       % Gy
        offaxis_idx = [1 3];   % Gx, Gz
    elseif contains(d, 'DIA')
        task_idx    = [1 2];   % Gx and Gy, each analyzed separately
        offaxis_idx = 3;       % Gz
    else
        % Unrecognized direction label -- don't assume an axis mapping.
        task_idx    = [1 2 3];
        offaxis_idx = [];
    end
end

function [peak_corr, peak_lag_sec] = axis_phase_relation(sig, ref, fs)
% Cross-correlates a head axis signal against the flange's direction-
% projected perturbation signal (rotation_native) to characterize timing.
% Returns the correlation value (sign preserved) at whichever lag has the
% largest |correlation|, and that lag converted to seconds. Requires
% xcorr() (Signal Processing Toolbox); falls back to NaN/NaN if it errors
% or the inputs are degenerate.
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
    max_lag = max(1, min(round(n/2), round(2 * fs)));  % search window: up to ~2s or half the trial
    try
        [c, lags] = xcorr(sig, ref, max_lag, 'coeff');
    catch
        return;  % e.g. xcorr not available
    end
    [~, idx] = max(abs(c));
    peak_corr    = c(idx);
    peak_lag_sec = lags(idx) / fs;
end

function [t_resamp, fs_resamp] = get_resampled_time_base(data)
% data.imu_*.six_axis_resamp is produced by the consumer via a normalized
% -length interp1 resample (Section 8b), not a direct decimation of the
% original clock, so there's no native per-sample timestamp to reuse.
% Reconstructs a synthetic time axis as sample-count / resample-rate,
% t = (1:N)/fs_resamp, using data.emg.fs (the rate the consumer actually
% resampled onto) when available, falling back to 2000 Hz (the driver's
% default emg_fs) only if that field is missing.
    N = size(data.imu_head.six_axis_resamp, 1);
    if isfield(data, 'emg') && isfield(data.emg, 'fs') && ~isempty(data.emg.fs) && isfinite(data.emg.fs)
        fs_resamp = data.emg.fs;
    else
        fs_resamp = 2000;  % fallback matching the driver's default emg_fs
    end
    t_resamp = (1:N)' / fs_resamp;
end

function m = compute_single_axis_metrics(sig, t, ref, fs)
% Computes magnitude (RMS, peak-to-peak, total excursion), directional
% bias (signed mean), and phase relation to the reference (flange
% perturbation) signal for one head gyro axis.
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

function block = compute_kinematics_block(head6, flange6, ref, fs, t, task_idx, offaxis_idx, suffix)
% Computes one full block of accel/gyro/off-axis-error metrics (mean/std
% for the 3 accel axes on both sensors, plus RMS/PtP/TotalExcursion/Bias/
% PhaseCorr/PhaseLag for each of the 3 head gyro axes, the off-axis-to-
% task RMS ratio, overall gyro-magnitude RMS, and duration), with every
% field name tagged by `suffix` (e.g. '' for native, '_Resamp' for the
% filtered/resampled data) so the two blocks can be merged into one row
% without colliding.
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

    head_gyro = head6(:, 4:6);   % [Gx Gy Gz]
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

    block.(['Head_Overall_GyroMag_RMS' suffix]) = sqrt(mean(sum(head_gyro.^2, 2), 'omitnan'));
    if numel(t) > 1
        block.(['Duration_sec' suffix]) = t(end) - t(1);
    else
        block.(['Duration_sec' suffix]) = NaN;
    end
end

function row = compute_trial_kinematics(data, tok, fpath)
% Computes flange/head per-axis average kinematics, plus the task-vs-
% off-axis compensation-error breakdown, for one trial -- TWICE: once on
% the native (Level-3, cycle-trimmed, native imu_fs) six-axis data, and
% once on the 5 Hz-lowpassed/resampled six-axis data (suffix '_Resamp'),
% via the same compute_kinematics_block() logic. Tagged with the input-
% parameter fields recovered from the filename.
    row = struct();

    % --- Linking / input-parameter fields ----------------------------------
    row.IMU_FileName = string(tok.imu_base);
    row.Muscle       = tok.muscle;
    row.Direction    = string(tok.direction);
    row.ParameterSet = string(tok.paramset);
    row.Condition    = string(tok.condition);
    row.Freq_Hz      = tok.freq_hz;
    row.Amp_raw      = tok.amp_raw;
    row.MatPath      = string(fpath);

    axis_names = {'Gx', 'Gy', 'Gz'};
    [task_idx, offaxis_idx] = axis_roles_for_direction(tok.direction);
    row.TaskAxes    = strjoin(axis_names(task_idx), ',');
    row.OffAxisAxes = strjoin(axis_names(offaxis_idx), ',');

    % --- Resampled block (5 Hz lowpass, resampled to EMG time base) --------
    % This is now the ONLY six-axis block used for head kinematics -- the
    % native block was retired along with data.imu_flange.six_axis_native,
    % which the consumer no longer saves (see its Section on space
    % savings). data.imu_head.six_axis_native itself is untouched and
    % still in the .mat file, but nothing here reads it anymore.
    [t_resamp, fs_resamp] = get_resampled_time_base(data);
    ref_resamp = get_flange_reference(data, tok.direction, 'six_axis_resamp');
    resamp_block = compute_kinematics_block(data.imu_head.six_axis_resamp, ...
        data.imu_flange.six_axis_resamp, ref_resamp, fs_resamp, t_resamp, ...
        task_idx, offaxis_idx, '_Resamp');

    row = merge_struct(row, resamp_block);

    % --- High-pass block (0.1 Hz zero-phase gyro high-pass, native imu_fs) -
    % Only present if this .mat file was produced by a consumer version
    % that saves six_axis_native_hp (accel unchanged, gyro high-pass
    % filtered on the full sync-trimmed signal before movement-buffer/
    % cycle-trim slicing -- see the consumer's Section 3b). Older .mat
    % files fall back to NaN-filled columns rather than erroring, so a
    % mixed-vintage batch still processes.
    %
    % NOTE on interpreting Head_G*_Bias_HP: a high-pass filter's entire
    % purpose is removing the near-DC component, so this Bias is expected
    % to land close to zero BY CONSTRUCTION -- it's a sanity check that
    % the filter worked, not an independent finding. The metrics that
    % actually matter here are RMS/PtP/TotalExcursion_HP -- the same
    % magnitude metrics as the (now-retired) native block, but with the
    % residual gyro calibration bias removed first.
    %
    % Reference: uses the flange's OWN high-pass-filtered six-axis data
    % (six_axis_native_hp) rather than the retired six_axis_native, so a
    % bias-corrected head signal is compared against an equivalently
    % bias-corrected flange reference rather than mixing in the retired
    % uncorrected one.
    if isfield(data.imu_head, 'six_axis_native_hp') && isfield(data.imu_flange, 'six_axis_native_hp')
        t_native  = data.imu_head.t_native(:);
        fs_native = data.imu_head.fs_native;
        ref_hp = get_flange_reference(data, tok.direction, 'six_axis_native_hp');
        hp_block = compute_kinematics_block(data.imu_head.six_axis_native_hp, ...
            data.imu_flange.six_axis_native_hp, ref_hp, fs_native, t_native, ...
            task_idx, offaxis_idx, '_HP');
        row.GyroHP_Available = true;
    else
        hp_block = nan_kinematics_block('_HP');
        row.GyroHP_Available = false;
    end
    row = merge_struct(row, hp_block);
end

function block = nan_kinematics_block(suffix)
% Same field set/names compute_kinematics_block would produce, all NaN --
% used when a .mat file predates the gyro high-pass fields, so the
% resulting table still has consistent columns across a mixed-vintage
% batch instead of erroring or silently dropping the block.
    accel_names = {'Ax', 'Ay', 'Az'};
    axis_names  = {'Gx', 'Gy', 'Gz'};
    block = struct();
    for c = 1:3
        block.(['Flange_' accel_names{c} '_mean' suffix]) = NaN;
        block.(['Flange_' accel_names{c} '_std' suffix])  = NaN;
    end
    for c = 1:3
        block.(['Flange_' axis_names{c} '_mean' suffix]) = NaN;
        block.(['Flange_' axis_names{c} '_std' suffix])  = NaN;
    end
    for c = 1:3
        block.(['Head_' accel_names{c} '_mean' suffix]) = NaN;
        block.(['Head_' accel_names{c} '_std' suffix])  = NaN;
    end
    for c = 1:3
        fn = axis_names{c};
        block.(['Head_' fn '_RMS' suffix])            = NaN;
        block.(['Head_' fn '_PtP' suffix])            = NaN;
        block.(['Head_' fn '_TotalExcursion' suffix]) = NaN;
        block.(['Head_' fn '_Bias' suffix])           = NaN;
        block.(['Head_' fn '_PhaseCorr' suffix])      = NaN;
        block.(['Head_' fn '_PhaseLag_sec' suffix])   = NaN;
    end
    block.(['OffAxis_to_Task_RMS_Ratio' suffix]) = NaN;
    block.(['Head_Overall_GyroMag_RMS' suffix])  = NaN;
    block.(['Duration_sec' suffix])              = NaN;
end

function s = merge_struct(s, extra)
% Copies every field of `extra` into `s` (overwriting on name collision).
    fn = fieldnames(extra);
    for i = 1:numel(fn)
        s.(fn{i}) = extra.(fn{i});
    end
end

function plot_trial_axis_vs_flange(t, head6, flange6, tok, out_dir, label)
% Saves a per-trial, 3-panel time-series figure -- one panel per axis
% (Gx, Gy, Gz) -- overlaying the head's signal on its matching flange
% axis, so the phase relationship PhaseCorr/PhaseLag summarize numerically
% can be inspected visually. `label` ('Native' or 'Resamp') selects the
% output subfolder and titles, so the same function serves both. Each
% panel is labeled TASK or OFF-AXIS per axis_roles_for_direction().
% Wrapped in try/catch (by the caller) so a plotting problem never breaks
% the batch loop.
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

    fname = sprintf('AxisVsFlange_%s_Muscle%d_%s_%s_%s_%s.png', ...
        tok.imu_base, tok.muscle, tok.direction, tok.paramset, tok.condition, label);
    out_png = fullfile(plot_dir, fname);

    try
        fig = figure('Visible', 'off', 'Position', [100 100 1000 800]);
        for c = 1:3
            subplot(3, 1, c);
            plot(t, flange6(:, 3 + c), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0, 'DisplayName', 'Flange');
            hold on;
            plot(t, head6(:, 3 + c), 'Color', [0.10 0.30 0.80], 'LineWidth', 1.0, 'DisplayName', 'Head');
            hold off;
            ylabel(sprintf('%s (native units)', axis_names{c}));
            title(sprintf('%s -- %s', axis_names{c}, role_label{c}));
            grid on;
            if c == 1
                legend('Location', 'best');
            end
            if c == 3
                xlabel('Time (s)');
            end
        end
        sgtitle(sprintf('%s | Muscle %d | %s | %s | %s | %s', ...
            tok.imu_base, tok.muscle, tok.direction, tok.paramset, tok.condition, label), ...
            'Interpreter', 'none');
        print(fig, out_png, '-dpng', '-r150');
        close(fig);
    catch ME
        fprintf('[WARNING] Could not generate %s axis-vs-flange plot for %s: %s\n', label, tok.imu_base, ME.message);
        if exist('fig', 'var') && isvalid(fig)
            close(fig);
        end
    end
end

function plot_head_rotation_overview(summary_table, out_dir)
% Saves diagnostic PNGs summarizing average head rotation and off-axis
% compensation error across the distinct ParameterSet/Condition/Direction
% combinations found. Uses the _Resamp columns -- the native (unsuffixed)
% block was retired along with data.imu_flange.six_axis_native. Purely
% diagnostic -- wrapped in try/catch so a plotting issue never breaks the
% pipeline.
    try
        labels = strcat(summary_table.ParameterSet, "_", ...
                         summary_table.Condition, "_", ...
                         summary_table.Direction);
        cat_labels = categorical(labels, labels);  % preserve original order

        % --- Fig 1: overall head rotation magnitude by condition ----------
        fig1 = figure('Visible', 'off', 'Position', [100 100 1200 500]);
        bar(cat_labels, summary_table.Head_Overall_GyroMag_RMS_Resamp);
        ylabel('Overall head gyro RMS (resampled, native units)');
        xlabel('ParameterSet\_Condition\_Direction');
        title('Average overall head rotation rate by parameter condition');
        xtickangle(45); grid on;
        out_png1 = fullfile(out_dir, 'Head_Rotation_Overview.png');
        print(fig1, out_png1, '-dpng', '-r150');
        close(fig1);

        % --- Fig 2: per-axis RMS (magnitude) by condition ------------------
        fig2 = figure('Visible', 'off', 'Position', [100 100 1200 500]);
        hold on;
        plot(cat_labels, summary_table.Head_Gx_RMS_Resamp, '-o', 'DisplayName', 'Gx RMS');
        plot(cat_labels, summary_table.Head_Gy_RMS_Resamp, '-o', 'DisplayName', 'Gy RMS');
        plot(cat_labels, summary_table.Head_Gz_RMS_Resamp, '-o', 'DisplayName', 'Gz RMS');
        hold off;
        legend('Location', 'best');
        ylabel('RMS (resampled, native units)');
        xlabel('ParameterSet\_Condition\_Direction');
        title('Average head gyro magnitude per axis by parameter condition');
        xtickangle(45); grid on;
        out_png2 = fullfile(out_dir, 'Head_Gyro_PerAxis_RMS.png');
        print(fig2, out_png2, '-dpng', '-r150');
        close(fig2);

        % --- Fig 3: off-axis compensation-error ratio by condition --------
        fig3 = figure('Visible', 'off', 'Position', [100 100 1200 500]);
        bar(cat_labels, summary_table.OffAxis_to_Task_RMS_Ratio_Resamp);
        ylabel('Off-axis / Task-axis RMS ratio (resampled)');
        xlabel('ParameterSet\_Condition\_Direction');
        title('Off-axis compensation error relative to task-axis movement');
        xtickangle(45); grid on;
        out_png3 = fullfile(out_dir, 'OffAxis_Error_Ratio.png');
        print(fig3, out_png3, '-dpng', '-r150');
        close(fig3);

        % --- Fig 4: directional bias per axis by condition -----------------
        fig4 = figure('Visible', 'off', 'Position', [100 100 1200 500]);
        bar(cat_labels, [summary_table.Head_Gx_Bias_Resamp, summary_table.Head_Gy_Bias_Resamp, summary_table.Head_Gz_Bias_Resamp]);
        legend({'Gx bias', 'Gy bias', 'Gz bias'}, 'Location', 'best');
        ylabel('Signed mean (resampled, native units)');
        xlabel('ParameterSet\_Condition\_Direction');
        title('Directional bias per axis by parameter condition');
        xtickangle(45); grid on; yline(0, 'k--');
        out_png4 = fullfile(out_dir, 'Head_Axis_Bias.png');
        print(fig4, out_png4, '-dpng', '-r150');
        close(fig4);

        fprintf('Diagnostic plots saved to:\n  %s\n  %s\n  %s\n  %s\n', ...
            out_png1, out_png2, out_png3, out_png4);
    catch ME
        fprintf('[WARNING] Could not generate summary plots: %s\n', ME.message);
    end
end