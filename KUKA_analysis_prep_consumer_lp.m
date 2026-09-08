%% ========================================================================
%  EMG / IMU ANALYSIS DATA PREPARATION
% ========================================================================
%  Builds on emg_imu_sync_pipeline_FIXED.m's sync + trim framework, but:
%
%   1) Processes ALL EMG channels (not just one verification channel) --
%      every channel gets: sync-trim -> rectify -> bandpass -> cycle-trim.
%
%   2) Explicitly organizes BOTH the flange IMU and the head IMU as
%      clearly labeled, parallel outputs. The head IMU is calibrated the
%      exact same way as the flange IMU (createCalibration + calibrateIMU
%      with its own static-pose calibration files), and trimmed the same
%      way: the sync-pulse trim (idxStart_imu:idxEnd_imu) applies to both
%      because they come from the same physical recording, and the
%      physical-swing cycle trim is estimated from the FLANGE signal only
%      (as in the original pipeline) and then applied identically to the
%      head IMU -- i.e. the head IMU's trim boundaries are DERIVED FROM,
%      and locked to, the flange IMU's cycle boundaries and the shared
%      sync signal, exactly like the current pipeline already does for
%      g2_lp. This script just makes that an explicit, top-level, named
%      output instead of an intermediate variable.
%
%   3) Packages everything (all EMG channels, both IMUs, resampled
%      signals, time bases, and processing metadata) into one struct
%      that's saved to disk, ready for downstream analysis scripts to
%      load without repeating any of the sync/trim/filter work.
%
%      EMG, per channel, several parallel processing levels (all derived
%      from the same >50 Hz high-pass signal):
%        data.emg.signal                    - high-pass -> rectify -> mu-centered bandpass -> cycle-trim (original "final" signal)
%        data.emg.signal_broadband_rect     - high-pass -> rectify -> cycle-trim (no mu-bandpass)
%        data.emg.signal_rms                - high-pass -> centered RMS (rms_win_sec, default 400ms) -> cycle-trim
%        data.emg.signal_highpass_cycletrim - high-pass -> cycle-trim only (no rectify/RMS/bandpass)
%        data.emg.signal_synctrim_highpass  - high-pass, sync-trim only (no cycle-trim)
%        data.emg.signal_synctrim_rect      - high-pass -> rectify, sync-trim only
%        data.emg.signal_synctrim_rms       - high-pass -> centered RMS, sync-trim only
%      IMU (flange/head/torso), in addition to the native + cycle-trimmed-
%      resampled fields: rotation_synctrim_resamp / six_axis_synctrim_resamp
%      / t_synctrim_resamp give the sync-trimmed (pre-cycle-trim) IMU
%      resampled onto the sync-trimmed EMG time base.
%
%      data.dprime.* (per EMG channel): a Gaussian-band-filtered (centered
%      at the swing frequency mu), forward-lag cross-correlation-aligned
%      (lag > 0 only, lag = 0 excluded), IMU-velocity-binned d-prime curve
%      -- see Section 8d for the full method, ported from dprime_2026.m.
%
%      All IMU lowpass filtering (the rotation_lp_cutoff_hz rotation-signal
%      filter, default 4 Hz, and the 6-axis imu_6axis_lp_cutoff_hz filter)
%      is applied ONCE to the FULL, untrimmed recording in Section 3, then
%      sliced at whatever trim boundaries each downstream section needs --
%      filtfilt never runs on an already-trimmed segment, so edge effects
%      don't bias samples near a trim/cycle boundary. Neither IMU rotation
%      signal is mean-subtracted anywhere in this file: MEMS calibration
%      (apply_mems_calibration, Section 3) already removes each gyro's DC
%      bias (b_gyro) before any of this runs, so re-subtracting a
%      per-trial mean downstream would double-correct the offset.
%
%  ------------------------------------------------------------------------
%  v2 -> v2 (functionized): this file used to be a script that assumed
%  rec_path, rec_num, imu_folder, calib1_file, calib2_file, imu_file,
%  emg_fs, imu_fs, emg_sync_ch, Muscle, direction, flange_calib_files,
%  head_calib_files, csv_out_dir already existed in the base workspace,
%  processing exactly one EMG-recording <-> IMU-file pair per run.
%
%  It is now a FUNCTION that takes all of those as fields of a single
%  input struct P, so it can be called in a loop -- once per row of the
%  sync_pulse_extractor's EMG_IMU_Match_Key.csv -- by a batch driver
%  script (see emg_imu_batch_data_prep.m). Calling it by hand for a
%  single pair still works exactly as before; just build P first:
%
%     P.rec_path  = rec_paths{1}; P.rec_num = 3; P.imu_folder = imu_folder;
%     P.calib1_file = calib1_file; P.calib2_file = calib2_file;
%     P.imu_file = 'KUKA_..._block1_2.txt'; P.emg_fs = emg_fs; ...
%     data = emg_imu_analysis_data_prep_v2(P);
%
%  Nothing about the processing itself (sync/trim/filter/calibration
%  logic below) has changed -- only how its inputs arrive and how a
%  sync-validation failure is reported (data.status instead of a bare
%  script `return`).
%
%  ------------------------------------------------------------------------
%  v2 -> v3 (sync indices now come from the extractor, not re-detected):
%  sync_extract_test_aug22.m already runs this exact hysteresis
%  low-high-low segmentation once per EMG recording / IMU file and saves
%  the resulting boundary sample indices in EMG_IMU_Matched_Trials.csv
%  (columns EMG_FirstHighToLowIdx / EMG_FinalLowToHighIdx /
%  IMU_FirstHighToLowIdx / IMU_FinalLowToHighIdx). This file used to
%  re-run that entire segmentation from scratch (duplicating
%  extract_all_segments/identify_trials as extract_sync_trial_hysteresis)
%  -- that's gone now. Instead, P must carry the four already-computed
%  indices for this EMG/IMU pair, straight from that CSV row:
%
%     P.idxStart_emg = row.EMG_FirstHighToLowIdx;
%     P.idxEnd_emg   = row.EMG_FinalLowToHighIdx;
%     P.idxStart_imu = row.IMU_FirstHighToLowIdx;
%     P.idxEnd_imu   = row.IMU_FinalLowToHighIdx;
%
%  The batch driver should join EMG_IMU_Matched_Trials.csv (on
%  EMG_SessionIdx/EMG_LocalRec/IMU_FileName, or simply iterate its rows
%  directly) to populate these before calling this function once per row.
%  P.emg_sync_ch is no longer used internally (kept optional, accepted
%  for backward compatibility with existing batch-driver P structs) --
%  the sync channel itself is never re-read here anymore.
% ========================================================================

function data = KUKA_analysis_prep_consumer_lp(P)
disp('>>> RUNNING NEW VERSION <<<');

rec_path           = P.rec_path;
rec_num            = P.rec_num;
imu_folder         = P.imu_folder;
calib1_file        = P.calib1_file;
calib2_file        = P.calib2_file;
imu_file           = P.imu_file;
emg_fs             = P.emg_fs;
imu_fs             = P.imu_fs;
Muscle             = P.Muscle;
direction          = P.direction;
flange_calib_files = P.flange_calib_files;
head_calib_files   = P.head_calib_files;
csv_out_dir        = P.csv_out_dir;
ParameterSet       = P.ParameterSet;   % from IMU_Key_Table, e.g. "1800_0504" -- used in the output filename AND as the swing-frequency estimate (see below)
condition          = P.condition;      % from IMU_Key_Table, "EO" or "EC" -- used only in the output filename

% The leading number in ParameterSet (e.g. "1320" in "1320_0300") is the
% swing-frequency estimate * 1000 -- i.e. divide by 1000 to get Hz
% (1320 -> 1.320 Hz). This replaces the old estimateNumCycles()-based
% frequency guess for cycle trimming and EMG bandpass centering (Section
% 6 below); NaN here means it couldn't be parsed and the code falls back
% to the old estimator with a warning.
freq_est_hz        = parse_freq_estimate_from_paramset(ParameterSet);

data = struct();
data.status = 'ok';   % overwritten below if sync validation fails
%% 0. IMU ROLE ASSIGNMENT (NEW -- explicit, top-level, parallel flange/head setup)
% imu_cal (built in Section 3 below) stacks each IMU's 6 calibrated
% columns [Ax Ay Az Gx Gy Gz] side by side in a fixed order:
%   columns  1:6  -> IMU 1 = FLANGE
%   columns  7:12 -> IMU 2 = HEAD
%   columns 13:18 -> IMU 3/4 = TORSO (only if present in this dataset)
% These index ranges are the single source of truth used everywhere below
% to slice out the flange vs. head IMU -- if your column order ever
% changes, this is the only place that needs to change with it.
FLANGE_IMU_COLS = 1:6;
HEAD_IMU_COLS   = 7:12;
TORSO_IMU_COLS  = 13:18;  % used only if a 3rd/4th IMU is detected

%% 1. LOAD DATA & EXTRACT STREAMS
session   = Session(rec_path);
node      = session.recordNodes{1};
recording = node.recordings{rec_num};
streams   = recording.continuous.keys();

if isempty(streams)
    error('No continuous streams found in recording %d.', rec_num);
end
stream = recording.continuous(streams{1});

tbl1     = load_imu_table_v2(fullfile(imu_folder, calib1_file));
tbl2     = load_imu_table_v2(fullfile(imu_folder, calib2_file));
imu_data = load_imu_table_v2(fullfile(imu_folder, imu_file));
fprintf(imu_folder)
% figure; plot(imu_data.Gx1); hold on; plot(imu_data.Gy1); plot(imu_data.Gz1);
% drawnow;

%% 2. SYNC PULSE INDICES (already computed by sync_extract_test_aug22.m)
% This used to re-run the extractor's whole hysteresis low-high-low
% segmentation from scratch (extract_sync_trial_hysteresis, ported from
% sync_pulse_extractor_v5_15th_pathfix_qwen_2.m). That's redundant: the
% extractor already ran it once per EMG recording / IMU file and saved
% the resulting boundary indices in EMG_IMU_Matched_Trials.csv
% (EMG_FirstHighToLowIdx/EMG_FinalLowToHighIdx/IMU_FirstHighToLowIdx/
% IMU_FinalLowToHighIdx -- exactly idxStart_emg/idxEnd_emg/idxStart_imu/
% idxEnd_imu below). Re-deriving them here risked silently diverging from
% the indices that were actually used to build the Match Key; instead we
% just take them as given.
if ~exist('dur_mismatch_warn_pct', 'var'); dur_mismatch_warn_pct = 2.0; end
if ~exist('dur_mismatch_fail_pct', 'var'); dur_mismatch_fail_pct = 10.0; end

required_idx_fields = {'idxStart_emg', 'idxEnd_emg', 'idxStart_imu', 'idxEnd_imu'};
missing_idx = required_idx_fields(~isfield(P, required_idx_fields));
if ~isempty(missing_idx)
    error(['P is missing required sync-trim index field(s): %s. These come ' ...
        'directly from EMG_IMU_Matched_Trials.csv (EMG_FirstHighToLowIdx / ' ...
        'EMG_FinalLowToHighIdx / IMU_FirstHighToLowIdx / IMU_FinalLowToHighIdx) ' ...
        '-- have the batch driver pass that row''s values through in P.'], ...
        strjoin(missing_idx, ', '));
end

idxStart_emg = P.idxStart_emg;   idxEnd_emg = P.idxEnd_emg;
idxStart_imu = P.idxStart_imu;   idxEnd_imu = P.idxEnd_imu;

emg_sync_ok = isfinite(idxStart_emg) && isfinite(idxEnd_emg) && idxEnd_emg > idxStart_emg;
imu_sync_ok = isfinite(idxStart_imu) && isfinite(idxEnd_imu) && idxEnd_imu > idxStart_imu;

if ~emg_sync_ok
    fprintf('[WARNING] EMG Sync: invalid/non-increasing indices supplied (idxStart=%g, idxEnd=%g).\n', ...
        idxStart_emg, idxEnd_emg);
end
if ~imu_sync_ok
    fprintf('[WARNING] IMU Sync: invalid/non-increasing indices supplied (idxStart=%g, idxEnd=%g).\n', ...
        idxStart_imu, idxEnd_imu);
end

dur_emg_sec  = (idxEnd_emg - idxStart_emg) / emg_fs;
dur_imu_sec  = (idxEnd_imu - idxStart_imu) / imu_fs;
dur_diff_pct = 100 * abs(dur_emg_sec - dur_imu_sec) / max(dur_imu_sec, eps);

fprintf('\n--- Initial Sync Trim Duration Check ---\n');
fprintf('EMG trimmed duration : %.4f s\n', dur_emg_sec);
fprintf('IMU trimmed duration : %.4f s\n', dur_imu_sec);
fprintf('Duration mismatch    : %.3f %%\n', dur_diff_pct);

if dur_diff_pct > dur_mismatch_fail_pct
    emg_sync_ok = false; imu_sync_ok = false;
    fprintf('[WARNING] Gross duration mismatch (%.2f%%). Treating sync as failed.\n', dur_diff_pct);
elseif dur_diff_pct > dur_mismatch_warn_pct
    fprintf('[WARNING] Duration mismatch (%.2f%%) exceeds soft tolerance -- proceeding, but verify.\n', dur_diff_pct);
end

if ~emg_sync_ok || ~imu_sync_ok
    fprintf('\n[WARNING] Trial skipped due to sync validation failure.\n');
    log_skipped_sync(csv_out_dir, imu_file, emg_sync_ok, imu_sync_ok);
    data.status = 'skipped_sync_fail';
    return;
end

%% 2b. MOVEMENT-START TRIM (LEVEL 2) -- fixed buffer inward from the Level-1
% sync-trim boundaries computed above. This is NOT the same thing as
% PreMarkerDur/PostMarkerDur -- those just describe the pre/post marker
% segments used to build the low-high-low triplet (and to match EMG
% trials to IMU trials in the Match Key). The movement-start buffer below
% is a fixed, separate amount of data trimmed off each end of the
% Level-1 sync trim to skip the participant's settle period before/after
% the actual movement. Override via P.movement_buffer_sec if needed.
if isfield(P, 'movement_buffer_sec')
    movement_buffer_sec = P.movement_buffer_sec;
else
    movement_buffer_sec = 4.0;
end
buf_imu = round(movement_buffer_sec * imu_fs);
buf_emg = round(movement_buffer_sec * emg_fs);

len_imu_sync = idxEnd_imu - idxStart_imu + 1;
len_emg_sync = idxEnd_emg - idxStart_emg + 1;
if len_imu_sync <= 2*buf_imu || len_emg_sync <= 2*buf_emg
    fprintf(['\n[WARNING] Movement buffer (%.2f s) leaves nothing after trimming ' ...
        '%.2f s off each end of a %.2f s (IMU) / %.2f s (EMG) sync-trimmed trial.\n'], ...
        movement_buffer_sec, movement_buffer_sec, len_imu_sync/imu_fs, len_emg_sync/emg_fs);
    log_skipped_sync(csv_out_dir, imu_file, emg_sync_ok, imu_sync_ok);
    data.status = 'skipped_movement_trim_fail';
    return;
end

idxStart_imu_mv = idxStart_imu + buf_imu;   idxEnd_imu_mv = idxEnd_imu - buf_imu;
idxStart_emg_mv = idxStart_emg + buf_emg;   idxEnd_emg_mv = idxEnd_emg - buf_emg;

% Relative offsets into the Level-1 (sync-trimmed) arrays built below --
% these are what's actually used to slice out the Level-2 dataset.
mv_start_rel_imu = buf_imu + 1;             mv_end_rel_imu = len_imu_sync - buf_imu;
mv_start_rel_emg = buf_emg + 1;             mv_end_rel_emg = len_emg_sync - buf_emg;

%% 3. IMU CALIBRATION (flange + head, identically processed) & INITIAL SYNC TRIM
% NOTE (ordering fix): the two-pose calibration vectors that DEFINE
% calibR_1/calibR_2/calibR_4 are now built from data that has been run
% through the SAME gyro-offset-removal + apply_mems_calibration
% correction the trial data gets (below, IMU1_AandG/IMU2_AandG), instead
% of raw accelerometer means. allori1/allori2 (tbl1/tbl2, the two static
% calibration poses) get their gyro DC offset removed the same way
% IMU1_AandG/IMU2_AandG do, before apply_mems_calibration runs on them.

try
    [imu1_c1_6ax, imu2_c1_6ax, imu4_c1_6ax] = extract_calib_arrays_6axis(tbl1);
catch ME
    error('Failed while parsing tbl1 (File: %s):\n%s', calib1_file, ME.message);
end
try
    [imu1_c2_6ax, imu2_c2_6ax, imu4_c2_6ax] = extract_calib_arrays_6axis(tbl2);
catch ME
    error('Failed while parsing tbl2 (File: %s):\n%s', calib2_file, ME.message);
end

% Each pose (allori1 = tbl1, allori2 = tbl2) is now run through the SAME
% gyro-offset removal + MEMS offset/scale/misalignment correction the
% trial data gets, using the matching per-IMU calib files, before its
% accelerometer mean is used to define the rotation matrix.
meanc1_imu1 = calibrated_pose_accel_mean(imu1_c1_6ax, flange_calib_files);
meanc2_imu1 = calibrated_pose_accel_mean(imu1_c2_6ax, flange_calib_files);
meanc1_imu2 = calibrated_pose_accel_mean(imu2_c1_6ax, head_calib_files);
meanc2_imu2 = calibrated_pose_accel_mean(imu2_c2_6ax, head_calib_files);

% Torso: there is currently no torso_calib_files field anywhere in the
% batch driver CONFIG. If a 3rd/4th IMU is actually in use, add one
% (P.torso_calib_files, same struct shape as flange/head) and pass it
% here -- and also apply it to IMU4_AandG in the trial-data section
% further down (it currently skips apply_mems_calibration entirely).
% Until then this falls back to the uncorrected mean, same as before.
if isfield(P, 'torso_calib_files')
    meanc1_imu4 = calibrated_pose_accel_mean(imu4_c1_6ax, P.torso_calib_files);
    meanc2_imu4 = calibrated_pose_accel_mean(imu4_c2_6ax, P.torso_calib_files);
else
    meanc1_imu4 = mean(imu4_c1_6ax(:, 1:3), 1, 'omitnan');
    meanc2_imu4 = mean(imu4_c2_6ax(:, 1:3), 1, 'omitnan');
end

calibR_1 = createCalibration(meanc1_imu1, meanc2_imu1);  % Flange
calibR_2 = createCalibration(meanc1_imu2, meanc2_imu2);  % Head
calibR_4 = createCalibration(meanc1_imu4, meanc2_imu4);  % Torso (if present)

vars = imu_data.Properties.VariableNames;

has_imu3_or_4 = any(ismember({'Ax4', 'Ax_4', 'Ax3', 'Ax_3'}, vars));
if ~has_imu3_or_4 && ~ismember('Ax1', vars) && ~ismember('Ax_1', vars)
    raw_arr = table2array(imu_data);
    if size(raw_arr, 2) >= 27
        has_imu3_or_4 = true;
    end
end
num_imus_detected = 2 + double(has_imu3_or_4);
fprintf('Detected %d IMU sensor(s) in current dataset.\n', num_imus_detected);

% --- FLANGE (IMU 1) ---
if ismember('Ax1', vars)
    IMU1_AandG = [imu_data.Ax1, imu_data.Ay1, imu_data.Az1, imu_data.Gx1, imu_data.Gy1, imu_data.Gz1];
elseif ismember('Ax_1', vars)
    IMU1_AandG = [imu_data.Ax_1, imu_data.Ay_1, imu_data.Az_1, imu_data.Gx_1, imu_data.Gy_1, imu_data.Gz_1];
else
    raw_arr    = table2array(imu_data);
    IMU1_AandG = raw_arr(:, 10:15);
end

% --- HEAD (IMU 2) ---
if ismember('Ax2', vars)
    IMU2_AandG = [imu_data.Ax2, imu_data.Ay2, imu_data.Az2, imu_data.Gx2, imu_data.Gy2, imu_data.Gz2];
elseif ismember('Ax_2', vars)
    IMU2_AandG = [imu_data.Ax_2, imu_data.Ay_2, imu_data.Az_2, imu_data.Gx_2, imu_data.Gy_2, imu_data.Gz_2];
else
    raw_arr    = table2array(imu_data);
    IMU2_AandG = raw_arr(:, 16:21);
end

% --- TORSO (IMU 3/4, if present) ---
if num_imus_detected == 3
    if ismember('Ax4', vars)
        IMU4_AandG = [imu_data.Ax4, imu_data.Ay4, imu_data.Az4, imu_data.Gx4, imu_data.Gy4, imu_data.Gz4];
    elseif ismember('Ax_4', vars)
        IMU4_AandG = [imu_data.Ax_4, imu_data.Ay_4, imu_data.Az_4, imu_data.Gx_4, imu_data.Gy_4, imu_data.Gz_4];
    elseif ismember('Ax3', vars)
        IMU4_AandG = [imu_data.Ax3, imu_data.Ay3, imu_data.Az3, imu_data.Gx3, imu_data.Gy3, imu_data.Gz3];
    elseif ismember('Ax_3', vars)
        IMU4_AandG = [imu_data.Ax_3, imu_data.Ay_3, imu_data.Az_3, imu_data.Gx_3, imu_data.Gy_3, imu_data.Gz_3];
    else
        raw_arr    = table2array(imu_data);
        IMU4_AandG = raw_arr(:, 22:27);
    end
end

% Gyro DC offset removal (identical treatment for every IMU)
% IMU1_AandG(:, 4:6) = IMU1_AandG(:, 4:6) - mean(IMU1_AandG(:, 4:6), 1, 'omitnan');
% IMU2_AandG(:, 4:6) = IMU2_AandG(:, 4:6) - mean(IMU2_AandG(:, 4:6), 1, 'omitnan');

fprintf('Applying MEMS matrix calibration to Flange IMU (IMU 1)...\n');
IMU1_AandG = apply_mems_calibration(IMU1_AandG, flange_calib_files);
fprintf('Applying MEMS matrix calibration to Head IMU (IMU 2)...\n');
IMU2_AandG = apply_mems_calibration(IMU2_AandG, head_calib_files);

imu1_cal = calibrateIMU(IMU1_AandG, calibR_1);  % Flange, calibrated
imu2_cal = calibrateIMU(IMU2_AandG, calibR_2);  % Head, calibrated

if num_imus_detected == 3
    % IMU4_AandG(:, 4:6) = IMU4_AandG(:, 4:6) - mean(IMU4_AandG(:, 4:6), 1, 'omitnan');
    imu4_cal = calibrateIMU(IMU4_AandG, calibR_4);
    imu_cal  = [imu1_cal, imu2_cal, imu4_cal];
else
    imu_cal  = [imu1_cal, imu2_cal];
end

% --- Full-recording lowpass filters (applied BEFORE any trimming) ---
% Both lowpass filters used downstream (the rotation_lp_cutoff_hz
% gyro-rotation lowpass, default 4 Hz, for direction projection/cycle
% detection in Section 5/6, and the 6-axis lowpass used ahead of every
% resampling/interpolation step in Sections 8b/8c) are built and applied here, to the FULL, untrimmed
% imu_cal, and only SLICED afterward at whatever indices each downstream
% section needs (sync trim, then cycle trim, or sync trim alone for the
% Section 8c sync-trim-only outputs). filtfilt is zero-phase but still
% has edge effects near the boundaries of whatever array it's given;
% running it once on the full recording and slicing afterward keeps
% every trim boundary's edge effect confined to the (much longer,
% discarded) untrimmed recording instead of biasing the samples right at
% the trim/cycle boundaries that downstream analysis actually uses.
if isfield(P, 'rotation_lp_cutoff_hz')
    rotation_lp_cutoff_hz = P.rotation_lp_cutoff_hz;
else
    rotation_lp_cutoff_hz = 4.0;
end
[b_lp, a_lp] = butter(2, rotation_lp_cutoff_hz / (imu_fs/2), 'low');   % rotation-signal LP, used for direction projection & cycle detection (Section 5/6)
if isfield(P, 'imu_6axis_lp_cutoff_hz')
    imu_6axis_lp_cutoff_hz = P.imu_6axis_lp_cutoff_hz;
else
    imu_6axis_lp_cutoff_hz = 5.0;
end
[b_lp6, a_lp6] = butter(2, imu_6axis_lp_cutoff_hz / (imu_fs/2), 'low');  % 6-axis LP, used ahead of six_axis_resamp / six_axis_synctrim_resamp (Section 8b/8c)

imu_cal_lp1 = filtfilt(b_lp, a_lp, imu_cal);    % rotation_lp_cutoff_hz LP, full recording, all IMU/axis columns
imu_cal_lp5 = filtfilt(b_lp6, a_lp6, imu_cal);  % imu_6axis_lp_cutoff_hz LP, full recording, all IMU/axis columns

% Initial sync trim -- identical for every IMU column, since they all
% come from the same physically-synced recording.
% idxStart_imu/idxEnd_imu now arrive pre-computed from
% EMG_IMU_Matched_Trials.csv rather than being derived from this file's
% own imu_cal, so bounds-check against the actual loaded length before
% slicing (a stale CSV row or mismatched file would otherwise throw an
% opaque indexing error here).
if idxStart_imu < 1 || idxEnd_imu > size(imu_cal, 1)
    error(['idxStart_imu:idxEnd_imu (%d:%d) from P falls outside the %d ' ...
        'samples loaded for %s -- the sync indices likely came from a ' ...
        'different IMU file/recording than P.imu_file.'], ...
        idxStart_imu, idxEnd_imu, size(imu_cal, 1), imu_file);
end
imu_cal_trimmed     = imu_cal(idxStart_imu:idxEnd_imu, :);
imu_cal_lp1_trimmed = imu_cal_lp1(idxStart_imu:idxEnd_imu, :);  % pre-filtered (full-recording), sync-trim slice
imu_cal_lp5_trimmed = imu_cal_lp5(idxStart_imu:idxEnd_imu, :);  % pre-filtered (full-recording), sync-trim slice
t_imu_trimmed   = (0:size(imu_cal_trimmed,1)-1) / imu_fs;

%% 4. MULTI-CHANNEL EMG EXTRACTION, SYNC TRIM & RECTIFICATION (ALL CHANNELS)
if Muscle == 1
    ch_range  = 1:64;
    ch_offset = 0;
elseif Muscle == 2
    ch_range  = 65:128;
    ch_offset = 64;
else
    error('Invalid Muscle index. Must be 1 or 2.');
end
n_emg_channels = length(ch_range);

raw_emg_array = double(stream.samples(ch_range, :));

[b_hp, a_hp] = butter(2, 50 / (emg_fs/2), 'high');
hp_emg_array = filtfilt(b_hp, a_hp, raw_emg_array')';

% Same bounds-check rationale as the IMU trim above -- idxStart_emg/
% idxEnd_emg came from the extractor's CSV, not from this recording.
if idxStart_emg < 1 || idxEnd_emg > size(hp_emg_array, 2)
    error(['idxStart_emg:idxEnd_emg (%d:%d) from P falls outside the %d ' ...
        'samples loaded for recording %d -- the sync indices likely came ' ...
        'from a different EMG recording than P.rec_path/P.rec_num.'], ...
        idxStart_emg, idxEnd_emg, size(hp_emg_array, 2), rec_num);
end
hp_emg_trimmed = hp_emg_array(:, idxStart_emg:idxEnd_emg);   % sync-trim, all channels
t_emg_trimmed  = (0:size(hp_emg_trimmed,2)-1) / emg_fs;
emg_rect_all_sync = abs(hp_emg_trimmed);                     % rectify, all channels

% --- Centered (non-causal) RMS envelope of the high-pass EMG ---
% Windowed RMS (default 400 ms) computed on the sync-trimmed, high-pass
% (>50 Hz) EMG -- i.e. sibling of emg_rect_all_sync above, but RMS instead
% of simple rectification. movmean's default centering (equal samples on
% each side of the current sample, 'Endpoints','shrink' to avoid NaNs/
% zero-padding bias at the very ends) makes this a centered, not causal,
% RMS. Override the window via P.emg_rms_window_sec if needed.
if isfield(P, 'emg_rms_window_sec')
    emg_rms_window_sec = P.emg_rms_window_sec;
else
    emg_rms_window_sec = 0.1;
end
rms_win_samples  = max(1, round(emg_rms_window_sec * emg_fs));
emg_rms_all_sync = sqrt(movmean(hp_emg_trimmed.^2, rms_win_samples, 2, 'Endpoints', 'shrink'));  % [n_emg_channels x N_sync], centered RMS

%% 5. FLANGE & HEAD ROTATION SIGNALS (direction handling, both IMUs in parallel)
% Raw (unfiltered), sync-trimmed gyro -- kept for the broadband signal
% further down (Section 5, "unfiltered counterparts").
gyro_flange = imu_cal_trimmed(:, FLANGE_IMU_COLS(4:6));
gyro_head   = imu_cal_trimmed(:, HEAD_IMU_COLS(4:6));

if num_imus_detected == 3
    gyro_torso = imu_cal_trimmed(:, TORSO_IMU_COLS(4:6));
else
    gyro_torso = zeros(size(gyro_flange));
end

% rotation_lp_cutoff_hz-lowpassed versions -- filtered on the FULL,
% untrimmed recording in Section 3 (imu_cal_lp1) and sliced here at the
% same sync-trim boundaries as gyro_flange/gyro_head/gyro_torso above, so
% filtfilt never runs on an already-trimmed (edge-sensitive) segment.
g_flange_lp_all = imu_cal_lp1_trimmed(:, FLANGE_IMU_COLS(4:6));
g_head_lp_all   = imu_cal_lp1_trimmed(:, HEAD_IMU_COLS(4:6));
if num_imus_detected == 3
    g_torso_lp_all = imu_cal_lp1_trimmed(:, TORSO_IMU_COLS(4:6));
else
    g_torso_lp_all = zeros(size(g_flange_lp_all));
end

% No mean subtraction here: the gyro DC offset was already removed by
% MEMS calibration (apply_mems_calibration subtracts calib_files.b_gyro,
% the per-IMU gyro bias, from every sample -- see Section 3 above and the
% apply_mems_calibration function below). Re-subtracting a per-trial mean
% on top of that would double-correct the offset and, because it's a
% windowed mean rather than the calibration's fixed bias, would also make
% the "zero" reference depend on how much of the trial is included in
% this particular array (sync-trim window here) rather than being a
% fixed, physically-anchored zero.
g_flange_xy = g_flange_lp_all(:, 1:2);
g_head_xy   = g_head_lp_all(:, 1:2);
g_torso_xy  = g_torso_lp_all(:, 1:2);

if strcmp(direction, 'ML')
    flange_rot = g_flange_xy(:, 1); head_rot = g_head_xy(:, 1); torso_rot = g_torso_xy(:, 1);
elseif strcmp(direction, 'AP')
    flange_rot = g_flange_xy(:, 2); head_rot = g_head_xy(:, 2); torso_rot = g_torso_xy(:, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %worked but offset by 1/4 cycle
% elseif contains(direction, 'Dia')
%     % PCA axis defined on the flange (the reference IMU), then the SAME
%     % rotation is applied to head/torso so all three signals are
%     % expressed in the same "swing direction" coordinate.
%     [coeff, ~, ~] = pca(g_flange_xy);
%     if coeff(2, 2) < 0
%         coeff(:, 2) = -coeff(:, 2);
%     end
%     flange_rot_full = g_flange_xy * coeff;
%     head_rot_full   = g_head_xy   * coeff;
%     torso_rot_full  = g_torso_xy  * coeff;
% 
%     flange_rot = flange_rot_full(:, 2);
%     head_rot   = head_rot_full(:, 2);
%     torso_rot  = torso_rot_full(:, 2);
elseif contains(direction, 'Dia')
    [coeff, ~, ~] = pca(g_flange_xy);          % PC1 = largest-variance axis

    % Choose the expected diagonal according to the trial label
    if contains(direction, 'R')                % DiaR, DiagonalR, etc.
        ref = [1;  1];                         % forward + right
    else                                       % DiaL, DiagonalL, …
        ref = [1; -1];                         % forward + left
    end
    ref = ref / norm(ref);

    % Flip PC1 so it points in the same half-plane as the reference
    if dot(coeff(:,1), ref) < 0
        coeff(:,1) = -coeff(:,1);
    end

    % (optional safety) if the two eigenvalues are almost equal,
    % you can fall back to the pure reference direction instead of PCA
    % [optional code omitted]

    flange_rot_full = g_flange_xy * coeff;
    head_rot_full   = g_head_xy   * coeff;
    torso_rot_full  = g_torso_xy  * coeff;

    flange_rot = flange_rot_full(:,1);         % signed PC1
    head_rot   = head_rot_full(:,1);
    torso_rot  = torso_rot_full(:,1);
else
    error('Unknown direction setting.');
end

Flange_LP_cut = flange_rot;   % kept for naming continuity with the sync/trim pipeline
Flange_LP_cut_full = Flange_LP_cut;

% --- Unfiltered, broadband counterparts, same projection ---
% flange_rot/head_rot above were built from g_flange_lp_all/g_head_lp_all,
% i.e. AFTER the rotation_lp_cutoff_hz lowpass -- appropriate for cycle
% detection & CRP, but it means any time-frequency plot of flange_rot/
% head_rot can never show content above that cutoff by construction.
% These broadband versions reuse the exact same axis/PCA projection, just
% applied to the RAW (calibrated, but unfiltered) gyro instead of the
% lowpassed gyro, so a time-frequency plot of them actually shows what
% frequencies are present. No mean subtraction here either, for the same
% reason as Section 5 above: apply_mems_calibration already removed the
% gyro's DC bias (b_gyro), so gyro_flange/gyro_head are already
% zero-referenced -- an additional per-trial mean subtraction would
% double-correct the offset.
gyro_flange_raw_xy = gyro_flange(:, 1:2);
gyro_head_raw_xy   = gyro_head(:, 1:2);

if strcmp(direction, 'ML')
    flange_rot_broadband = gyro_flange_raw_xy(:, 1);
    head_rot_broadband   = gyro_head_raw_xy(:, 1);
elseif strcmp(direction, 'AP')
    flange_rot_broadband = gyro_flange_raw_xy(:, 2);
    head_rot_broadband   = gyro_head_raw_xy(:, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%works but like other section, it offsets by 1/4 cycle
% elseif contains(direction, 'Dia')
%     flange_rot_broadband_full = gyro_flange_raw_xy * coeff;   % same coeff computed above from the LOWPASSED flange signal, just applied to raw data here
%     head_rot_broadband_full   = gyro_head_raw_xy   * coeff;
%     flange_rot_broadband = flange_rot_broadband_full(:, 2);
%     head_rot_broadband   = head_rot_broadband_full(:, 2);
elseif contains(direction, 'Dia')
    flange_rot_broadband_full = gyro_flange_raw_xy * coeff;   % same coeff
    head_rot_broadband_full   = gyro_head_raw_xy   * coeff;
    flange_rot_broadband = flange_rot_broadband_full(:, 1);   % PC1
    head_rot_broadband   = head_rot_broadband_full(:, 1);
end

%% 6. PHYSICAL SWING ESTIMATION & CYCLE TRIM (estimated from FLANGE, applied to everything)
duration_imu_sec = length(Flange_LP_cut) / imu_fs;
num_cycles       = estimateNumCycles(Flange_LP_cut);   % diagnostic/reporting only -- see note below

% Prefer the frequency estimate carried in ParameterSet (IMU_Key_Table)
% over estimateNumCycles(). estimateNumCycles() counts swings across the
% WHOLE post-sync signal, including the movement_buffer_sec settle period
% at each end where little/no swinging happens; dividing by the full
% duration therefore underestimates frequency (overestimates the cycle
% period), which can make Section 6's one-cycle trim larger than the
% active window itself. mu/T_cycle here also set the EMG movement-locked
% bandpass center further down, so getting this right matters beyond
% just the cycle trim.
if ~isnan(freq_est_hz) && freq_est_hz > 0
    mu      = freq_est_hz;
    T_cycle = 1 / mu;
    fprintf('\nUsing ParameterSet frequency estimate: %.4f Hz (Period T = %.4f s) [from "%s"]\n', ...
        mu, T_cycle, ParameterSet);
else
    mu      = num_cycles / duration_imu_sec;
    T_cycle = 1 / mu;
    warning(['ParameterSet ("%s") did not yield a usable frequency estimate; falling back to ' ...
        'the swing-count-over-full-duration estimate (%.4f Hz). This estimate is unreliable ' ...
        'when swings are concentrated in the active (non-buffer) window.'], ...
        char(ParameterSet), mu);
end

fprintf('\n--- Physical Swing & Cycle Trim Diagnostic ---\n');
fprintf('Trimmed Duration (Post-Sync) : %.2f seconds\n', duration_imu_sec);
fprintf('Detected Physical Swings     : %.1f swings (%.2f cycles) [diagnostic only]\n', num_cycles * 2, num_cycles);
fprintf('Physical Frequency Used      : %.4f Hz (Period T = %.4f s)\n', mu, T_cycle);

nom_trim_imu      = round(T_cycle * imu_fs);
debounce_samples  = max(1, round(0.15 * nom_trim_imu));

fprintf('Movement buffer : %.2f s (%d samples) trimmed off each end before cycle search.\n', ...
    movement_buffer_sec, buf_imu);

% --- Cycle trim via ACTUAL zero-crossings on the main (flange) axis ---
% Crossings alternate direction (rising, falling, rising, ...), so any
% two crossings 2 positions apart are the same direction and therefore
% ALWAYS exactly one period apart, regardless of where in the cycle the
% first one happens to fall. That's guaranteed -- but where the first
% crossing (c1) itself lands relative to the window edge is NOT fixed:
% it depends on the signal's phase when the movement-buffer window opens.
%   - If the signal is near zero right at the window edge, c1 falls
%     close to the edge (near 0 samples in) -- c1 IS close to a
%     "phase-zero" start, so counting 2 crossings from it (c1 -> c3, one
%     full period) lands close to exactly one period from the edge.
%   - If the signal is already well off zero at the window edge (e.g.
%     just above zero and swinging up first, as in the case that
%     motivated this), the search has to wait for the signal to swing up,
%     over the peak, and back down before finding ANY crossing -- so c1
%     itself falls near T/2 (half a period) in. In that case, c1 is
%     really "the half-cycle mark," and counting 2 crossings from it
%     (c1 -> c3) overshoots to 1.5 periods from the edge; counting only 1
%     crossing (c1 -> c2) is what actually lands near one period.
%
% So which count (1 or 2 crossings past c1) is correct depends on where
% c1 itself falls -- and the natural threshold for that decision is
% T/4: if c1 is within a quarter-period of the window edge, treat it as
% the phase-zero case (use c3); otherwise treat it as the half-cycle-mark
% case (use c2). This is applied symmetrically at the end of the window
% using the LAST detected crossing and its distance from the end edge.
all_crossings = find_debounced_zero_cross(Flange_LP_cut, mv_start_rel_imu, mv_end_rel_imu, debounce_samples, 'all');

fprintf('Zero-crossings found in active window : %d\n', numel(all_crossings));

quarter_period_imu = nom_trim_imu / 4;   % T/4 in samples, from the ParameterSet/estimateNumCycles frequency

if numel(all_crossings) >= 5
    % --- Start boundary ---
    elapsed_to_c1 = all_crossings(1) - mv_start_rel_imu;
    if elapsed_to_c1 < quarter_period_imu
        idxStart_imu_cycle = all_crossings(3);   % c1 ~ phase-zero -> one full period is c1->c3
        start_rule = 'c1 within T/4 of window start -> used c3 (c1 + 1 full cycle)';
    else
        idxStart_imu_cycle = all_crossings(2);   % c1 ~ half-cycle mark -> one full period is c1->c2
        start_rule = 'c1 beyond T/4 of window start -> used c2 (c1 + half cycle)';
    end

    % --- End boundary (mirrored) ---
    elapsed_from_last = mv_end_rel_imu - all_crossings(end);
    if elapsed_from_last < quarter_period_imu
        idxEnd_imu_cycle = all_crossings(end-2);
        end_rule = 'last crossing within T/4 of window end -> used end-2 (last - 1 full cycle)';
    else
        idxEnd_imu_cycle = all_crossings(end-1);
        end_rule = 'last crossing beyond T/4 of window end -> used end-1 (last - half cycle)';
    end

    fprintf('Cycle-trim start: %s (elapsed=%d samples, T/4=%.1f samples)\n', start_rule, elapsed_to_c1, quarter_period_imu);
    fprintf('Cycle-trim end  : %s (elapsed=%d samples, T/4=%.1f samples)\n', end_rule, elapsed_from_last, quarter_period_imu);

    if idxEnd_imu_cycle <= idxStart_imu_cycle
        % Degenerate case: only possible with a very low crossing count
        % (e.g. numel==5, where both the "c3" and "end-2" choices can
        % land on the same middle crossing). Fall back rather than error.
        idxStart_imu_cycle = mv_start_rel_imu + nom_trim_imu;
        idxEnd_imu_cycle   = mv_end_rel_imu   - nom_trim_imu;
        warning(['Cycle-trim start/end landed on the same crossing (or crossed) -- too few ' ...
            'crossings (%d found) for the counting rule to resolve both boundaries ' ...
            'independently. Falling back to a fixed one-period offset (%.4f Hz, %d samples) ' ...
            'from the movement-buffer boundaries.'], numel(all_crossings), mu, nom_trim_imu);
    end
else
    % Not enough real crossings to apply the counting rule on each side --
    % fall back to the fixed one-period offset (still using the
    % ParameterSet-derived T_cycle, so at least the period itself is
    % trustworthy even though its placement isn't verified against the
    % actual signal).
    idxStart_imu_cycle = mv_start_rel_imu + nom_trim_imu;
    idxEnd_imu_cycle   = mv_end_rel_imu   - nom_trim_imu;
    warning(['Only %d zero-crossing(s) found in the active window -- not enough to apply the ' ...
        'T/4 counting rule on each side (need >=5). Falling back to a fixed one-period offset ' ...
        '(%.4f Hz, %d samples) from the movement-buffer boundaries.'], ...
        numel(all_crossings), mu, nom_trim_imu);
end

if idxEnd_imu_cycle <= idxStart_imu_cycle
    error('Cycle trim failed: idxEnd_imu_cycle (%d) <= idxStart_imu_cycle (%d).', ...
        idxEnd_imu_cycle, idxStart_imu_cycle);
end

% Crop the flange, head, and torso rotation signals at the SAME
% flange-derived cycle boundaries -- this is the "head IMU trimmed
% identically to the flange IMU" step.
flange_rot_cycle = flange_rot(idxStart_imu_cycle : idxEnd_imu_cycle);
head_rot_cycle   = head_rot(idxStart_imu_cycle : idxEnd_imu_cycle);
torso_rot_cycle  = torso_rot(idxStart_imu_cycle : idxEnd_imu_cycle);
Flange_LP_cut    = flange_rot_cycle(:);

flange_rot_broadband_cycle = flange_rot_broadband(idxStart_imu_cycle : idxEnd_imu_cycle);
head_rot_broadband_cycle   = head_rot_broadband(idxStart_imu_cycle : idxEnd_imu_cycle);

% Also crop the FULL 6-axis calibrated flange/head (and torso) streams at
% the same boundaries, so downstream analyses have access to more than
% just the single swing-direction rotation signal if needed.
imu_flange_6axis_cycle = imu_cal_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, FLANGE_IMU_COLS);
imu_head_6axis_cycle   = imu_cal_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, HEAD_IMU_COLS);
if num_imus_detected == 3
    imu_torso_6axis_cycle = imu_cal_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, TORSO_IMU_COLS);
end
t_imu_cycle = (0:length(Flange_LP_cut)-1) / imu_fs;

% Map cycle-trim boundaries onto the EMG time base (ratio mapping, same
% approach as the sync/trim pipeline)
ratio_start = (idxStart_imu_cycle - 1) / (length(imu_cal_trimmed) - 1);
ratio_end   = (idxEnd_imu_cycle - 1)   / (length(imu_cal_trimmed) - 1);

%% 7. BANDPASS ALL EMG CHANNELS, THEN CYCLE-TRIM (sync-trim -> rectify -> bandpass -> cycle-trim)
half_bw = 0.05;
f_low   = max(0.03, mu - half_bw);
f_high  = mu + half_bw;
[b_bp, a_bp] = butter(1, [f_low, f_high] / (emg_fs/2), 'bandpass');

emg_bp_all_sync = filtfilt(b_bp, a_bp, emg_rect_all_sync')';   % ALL channels, full sync-trim length

% --- Sync-trim-ONLY versions (no cycle trim), for diagnostic comparison ---
% emg_bp_all_sync above is already exactly this for EMG: sync-trimmed,
% rectified, mu-bandpassed, but NOT yet cropped to the physical-swing
% cycle boundaries. Here we resample the (also not-yet-cycle-trimmed)
% flange_rot/head_rot onto that same EMG time base, so you can compute
% CRP etc. over the FULL sync-trimmed trial, including whatever the
% cycle-trim step in Section 6 cuts off the front/back.
N_sync_only = size(emg_bp_all_sync, 2);
Flange_rot_sync_only_resamp = interp1(linspace(0,1,length(flange_rot)), flange_rot, linspace(0,1,N_sync_only), 'linear')';
Head_rot_sync_only_resamp   = interp1(linspace(0,1,length(head_rot)),   head_rot,   linspace(0,1,N_sync_only), 'linear')';
t_sync_only = (0:N_sync_only-1) / emg_fs;

idxStart_emg_cycle = 1 + round(ratio_start * (size(emg_bp_all_sync, 2) - 1));
idxEnd_emg_cycle   = 1 + round(ratio_end   * (size(emg_bp_all_sync, 2) - 1));

EMG_all = emg_bp_all_sync(:, idxStart_emg_cycle : idxEnd_emg_cycle);  % ALL channels, final
t_emg_final = (0:size(EMG_all,2)-1) / emg_fs;

% Broadband EMG: high-pass (>50 Hz) + rectified only, i.e. WITHOUT the
% mu +/- 0.1 Hz movement-locked bandpass above, cycle-trimmed to the same
% boundaries as EMG_all. The narrow mu-centered bandpass in EMG_all is
% appropriate for amplitude-modulation analyses like continuous relative
% phase, but it removes essentially all the 20-450 Hz myoelectric content
% a genuine EMG time-frequency (spectrogram/scalogram) plot needs. Keep
% this broadband version around for that purpose.
EMG_broadband_rect_all = emg_rect_all_sync(:, idxStart_emg_cycle : idxEnd_emg_cycle);  % ALL channels

% Cycle-trimmed RMS envelope: same idxStart_emg_cycle:idxEnd_emg_cycle
% boundaries as EMG_broadband_rect_all above, applied to the centered-RMS
% array computed in Section 4 -- sync-trim -> high-pass -> RMS(400ms) ->
% cycle-trim, i.e. the RMS sibling of EMG_broadband_rect_all/EMG_all.
EMG_rms_all = emg_rms_all_sync(:, idxStart_emg_cycle : idxEnd_emg_cycle);  % ALL channels, centered RMS, cycle-trimmed

% Cycle-trimmed, high-pass-ONLY EMG (no rectify, no mu-centered bandpass):
% sync-trim -> high-pass -> cycle-trim, same boundaries as above.
EMG_highpass_cycle_all = hp_emg_trimmed(:, idxStart_emg_cycle : idxEnd_emg_cycle);  % ALL channels

dur_imu_cycle_sec = (idxEnd_imu_cycle - idxStart_imu_cycle) / imu_fs;
dur_emg_cycle_sec = (idxEnd_emg_cycle - idxStart_emg_cycle) / emg_fs;
dur_cycle_diff_pct = 100 * abs(dur_imu_cycle_sec - dur_emg_cycle_sec) / max(dur_imu_cycle_sec, eps);
fprintf('\n--- Cycle Trim Duration Check ---\n');
fprintf('IMU cycle-trimmed duration : %.4f s\n', dur_imu_cycle_sec);
fprintf('EMG cycle-trimmed duration : %.4f s\n', dur_emg_cycle_sec);
fprintf('Duration mismatch          : %.3f %%\n', dur_cycle_diff_pct);
if dur_cycle_diff_pct > dur_mismatch_warn_pct
    fprintf('[WARNING] Cycle-trim duration mismatch exceeds %.1f%%.\n', dur_mismatch_warn_pct);
end

%% 8. RESAMPLE FLANGE & HEAD ROTATION SIGNALS TO THE FINAL EMG TIME BASE
N_final = size(EMG_all, 2);
Flange_rot_2000Hz = interp1(linspace(0,1,length(flange_rot_cycle)), flange_rot_cycle, linspace(0,1,N_final), 'linear')';
Head_rot_2000Hz   = interp1(linspace(0,1,length(head_rot_cycle)),   head_rot_cycle,   linspace(0,1,N_final), 'linear')';
if num_imus_detected == 3
    Torso_rot_2000Hz = interp1(linspace(0,1,length(torso_rot_cycle)), torso_rot_cycle, linspace(0,1,N_final), 'linear')';
end

%% 8b. FULL 6-AXIS IMU: RESAMPLE THE (ALREADY FULL-RECORDING-LOWPASSED)
% SIGNAL ONTO THE SAME FINAL TIME BASE AS THE ROTATION SIGNALS ABOVE
% (i.e. matches EMG_all's N_final samples @ emg_fs). Each of the 6
% columns (Ax Ay Az Gx Gy Gz) was already run through a bidirectional
% (filtfilt, zero-phase) 2nd-order Butterworth lowpass at
% imu_6axis_lp_cutoff_hz (default 5 Hz) BACK IN SECTION 3, on the FULL,
% untrimmed recording (imu_cal_lp5) -- same filter family/order as the
% rotation_lp_cutoff_hz gyro lowpass used for cycle detection, just a
% different cutoff.
% Here we only SLICE that pre-filtered array at the cycle-trim boundaries
% (idxStart_imu_cycle:idxEnd_imu_cycle, relative to the sync-trimmed
% array, same indexing as imu_cal_lp5_trimmed / imu_flange_6axis_cycle
% etc.) -- filtfilt is never run on an already-trimmed segment, so trim
% boundaries don't bias the filtered samples nearest to them. Resampling
% uses the same normalized-length interp1 approach used for
% flange_rot_cycle/head_rot_cycle just above (NOT a raw-time
% t_imu_cycle->t_emg_final mapping), so a small native duration mismatch
% between the IMU and EMG cycle-trimmed lengths doesn't require
% extrapolation at the edges.
imu_flange_6axis_lp = imu_cal_lp5_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, FLANGE_IMU_COLS);
imu_head_6axis_lp   = imu_cal_lp5_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, HEAD_IMU_COLS);

resamp_6axis = @(sig6) cell2mat(arrayfun(@(c) ...
    interp1(linspace(0,1,size(sig6,1)), sig6(:,c), linspace(0,1,N_final), 'linear')', ...
    1:size(sig6,2), 'UniformOutput', false));

Flange_6axis_2000Hz = resamp_6axis(imu_flange_6axis_lp);   % [N_final x 6]
Head_6axis_2000Hz   = resamp_6axis(imu_head_6axis_lp);     % [N_final x 6]

if num_imus_detected == 3
    imu_torso_6axis_lp = imu_cal_lp5_trimmed(idxStart_imu_cycle:idxEnd_imu_cycle, TORSO_IMU_COLS);
    Torso_6axis_2000Hz = resamp_6axis(imu_torso_6axis_lp);  % [N_final x 6]
end

%% 8c. SYNC-TRIMMED (NOT cycle-trimmed) IMU, RESAMPLED ONTO THE SYNC-TRIMMED
% EMG TIME BASE. Same idea as Section 8/8b, but for the Level-1 sync-trim
% length (N_sync_only samples, i.e. matching emg_bp_all_sync /
% data.emg.signal_synctrim_* before any cycle trim) instead of the final
% cycle-trimmed length (N_final). Flange_rot_sync_only_resamp /
% Head_rot_sync_only_resamp were already computed just above (Section 7)
% for CRP-style use; this just adds the matching torso rotation and the
% full 6-axis blocks (same 5 Hz lowpass as imu_6axis_lp_cutoff_hz, applied
% at native imu_fs before resampling) for completeness.
if num_imus_detected == 3
    Torso_rot_sync_only_resamp = interp1(linspace(0,1,length(torso_rot)), torso_rot, linspace(0,1,N_sync_only), 'linear')';
end

resamp_6axis_sync = @(sig6) cell2mat(arrayfun(@(c) ...
    interp1(linspace(0,1,size(sig6,1)), sig6(:,c), linspace(0,1,N_sync_only), 'linear')', ...
    1:size(sig6,2), 'UniformOutput', false));

% Same fix as Section 8b: slice the FULL-RECORDING-lowpassed imu_cal_lp5
% (built once in Section 3, before any trimming) at the sync-trim
% boundaries, rather than filtering the already sync-trimmed
% imu_cal_trimmed segment -- keeps filtfilt off of trimmed/edge-sensitive
% data here too. imu_cal_lp5_trimmed is exactly this: imu_cal_lp5 sliced
% at idxStart_imu:idxEnd_imu, i.e. the sync-trim length.
imu_flange_6axis_synctrim_lp = imu_cal_lp5_trimmed(:, FLANGE_IMU_COLS);
imu_head_6axis_synctrim_lp   = imu_cal_lp5_trimmed(:, HEAD_IMU_COLS);

Flange_6axis_synctrim_resamp = resamp_6axis_sync(imu_flange_6axis_synctrim_lp);   % [N_sync_only x 6]
Head_6axis_synctrim_resamp   = resamp_6axis_sync(imu_head_6axis_synctrim_lp);     % [N_sync_only x 6]

if num_imus_detected == 3
    imu_torso_6axis_synctrim_lp = imu_cal_lp5_trimmed(:, TORSO_IMU_COLS);
    Torso_6axis_synctrim_resamp = resamp_6axis_sync(imu_torso_6axis_synctrim_lp);  % [N_sync_only x 6]
end

%% 8d. GAUSSIAN-FILTERED CROSS-CORRELATION & D-PRIME (PER EMG CHANNEL)
% Ported from dprime_2026.m (Section 5: "PER-CHANNEL GAUSSIAN FILTERING,
% CROSS-CORRELATION, & D-PRIME"). For each EMG channel independently:
%   1) A narrow Gaussian band-filter, centered at +/-mu (the swing
%      frequency already established in Section 6), is built in the
%      frequency domain and applied (via FFT multiply / IFFT) to both the
%      rectified EMG amplitude (EMG_broadband_rect_all) and the flange
%      velocity signal (Flange_rot_2000Hz) -- both already share the same
%      N_final-sample, emg_fs time base, so no separate resampling is
%      needed here.
%   2) The channel's filtered EMG is cross-correlated against the
%      filtered IMU to find that channel's own best lag, and the two
%      signals are aligned at that lag.
%   3) The aligned EMG is binned by the aligned IMU velocity
%      (bin width dprime_bin_width, in deg/s), and d-prime is computed
%      between each velocity bin and the near-zero-velocity "off" bin:
%         d' = |mean_on - mean_off| / sqrt((var_on + var_off) / 2)
%
% CROSS-CORRELATION SEARCH -- forward-looking, lag=0 excluded:
% Only STRICTLY POSITIVE lags are searched (0 < lag < dprime_max_lag_sec),
% i.e. only alignments where the Gaussian-filtered EMG must be shifted
% forward in time to line up with the IMU are considered. This differs
% from the original dprime_2026.m, which allowed lag = 0; excluding it
% here prevents the search from trivially "aligning" at no shift and
% forces it to find a genuine (nonzero) lag. Override the search window
% via P.dprime_xcorr_max_lag_sec (default 2 s, same as the original
% script's two_sec_limit).
%
% All parameters are overridable via P (all optional, defaults shown):
%   P.dprime_fwhm_hz              (0.01)  Gaussian filter FWHM, Hz
%   P.dprime_bin_width_deg_s      (0.20)  velocity bin width, deg/s
%   P.dprime_xcorr_max_lag_sec    (2.0)   max forward lag searched, s
%   P.dprime_min_samples_per_bin  (20)    min samples for a bin to be used

if isfield(P, 'dprime_fwhm_hz')
    dprime_fwhm = P.dprime_fwhm_hz;
else
    dprime_fwhm = 0.01;
end
if isfield(P, 'dprime_bin_width_deg_s')
    dprime_bin_width = P.dprime_bin_width_deg_s;
else
    dprime_bin_width = 0.20;
end
if isfield(P, 'dprime_xcorr_max_lag_sec')
    dprime_max_lag_sec = P.dprime_xcorr_max_lag_sec;
else
    dprime_max_lag_sec = 2.0;
end
if isfield(P, 'dprime_min_samples_per_bin')
    dprime_min_samples_per_bin = P.dprime_min_samples_per_bin;
else
    dprime_min_samples_per_bin = 20;
end

% --- Gaussian filter, built at N_final/emg_fs resolution, centered at mu ---
dprime_sigma = dprime_fwhm / (2 * sqrt(2 * log(2)));
df_dprime    = emg_fs / N_final;
f_dprime     = (0:N_final-1) * df_dprime;
f2_dprime    = f_dprime - emg_fs * (f_dprime >= emg_fs/2);

gaussian_filter_pos    = exp(-((f2_dprime - mu).^2) / (2 * dprime_sigma^2))';
gaussian_filter_neg    = exp(-((f2_dprime - (-mu)).^2) / (2 * dprime_sigma^2))';
dprime_gaussian_filter  = gaussian_filter_pos + gaussian_filter_neg;
dprime_gaussian_filter  = dprime_gaussian_filter / max(dprime_gaussian_filter);

% --- Apply filter to the (shared) IMU velocity signal once ---
imu_fft_dprime       = fft(Flange_rot_2000Hz);
filtered_IMU_dprime  = ifft(imu_fft_dprime .* dprime_gaussian_filter, 'symmetric'); %#ok<NASGU> % kept for diagnostics/plots

% --- Velocity bin setup (shared across channels) ---
dprime_max_val          = max(abs(Flange_rot_2000Hz));
dprime_num_bins_per_side = ceil(dprime_max_val / dprime_bin_width);
dprime_bin_centers  = (-dprime_num_bins_per_side : dprime_num_bins_per_side) * dprime_bin_width;
dprime_bin_edges    = dprime_bin_centers - dprime_bin_width/2;
dprime_bin_edges    = [dprime_bin_edges, dprime_bin_centers(end) + dprime_bin_width/2];
num_dprime_bins     = length(dprime_bin_centers);
dprime_off_bin_idx  = find(abs(dprime_bin_centers) < 1e-5, 1);

dprime_all      = nan(n_emg_channels, num_dprime_bins);
dprime_lags_sec = nan(n_emg_channels, 1);

if isempty(dprime_off_bin_idx)
    warning(['dprime: no velocity bin centered at 0 deg/s -- there is no baseline ' ...
        '"off" bin, so d-prime cannot be computed for any channel. Leaving ' ...
        'data.dprime.values as all-NaN.']);
else
    dprime_xcorr_max_lag_samp = round(dprime_max_lag_sec * emg_fs);

    for ch = 1:n_emg_channels
        emg_ch          = EMG_broadband_rect_all(ch, :)';
        emg_fft_ch       = fft(emg_ch);
        filtered_EMG_ch  = ifft(emg_fft_ch .* dprime_gaussian_filter, 'symmetric');

        % Forward-looking cross-correlation: search strictly positive lags
        % only (lag = 0 excluded), within +/- dprime_max_lag_sec.
        [xc, xc_lags] = xcorr(Flange_rot_2000Hz, filtered_EMG_ch);
        valid_lags = (xc_lags > 0) & (xc_lags < dprime_xcorr_max_lag_samp);
        if ~any(valid_lags)
            warning(['dprime: channel %d (global ch %d) has no valid strictly-positive ' ...
                'lag within %.2f s -- skipping this channel (left as NaN).'], ...
                ch, ch_range(ch), dprime_max_lag_sec);
            continue;
        end
        pos_corr = xc(valid_lags);
        pos_lags = xc_lags(valid_lags);

        [~, maxidx]  = max(abs(pos_corr));
        optimal_lag  = pos_lags(maxidx);
        dprime_lags_sec(ch) = optimal_lag / emg_fs;

        emg_aligned = filtered_EMG_ch(1:end-optimal_lag);
        imu_aligned = Flange_rot_2000Hz(optimal_lag+1:end);

        % Velocity binning for this channel
        bin_mean = nan(1, num_dprime_bins);
        bin_var  = nan(1, num_dprime_bins);
        for b = 1:num_dprime_bins
            in_bin = find(imu_aligned >= dprime_bin_edges(b) & imu_aligned < dprime_bin_edges(b+1));
            if length(in_bin) >= dprime_min_samples_per_bin
                bin_mean(b) = mean(emg_aligned(in_bin));
                bin_var(b)  = var(emg_aligned(in_bin));
            end
        end

        mu_off  = bin_mean(dprime_off_bin_idx);
        var_off = bin_var(dprime_off_bin_idx);
        for b = 1:num_dprime_bins
            if isnan(bin_mean(b)) || isnan(mu_off)
                dprime_all(ch, b) = NaN;
            else
                dprime_all(ch, b) = abs(bin_mean(b) - mu_off) / sqrt((bin_var(b) + var_off) / 2);
            end
        end
    end
end

%% 9. PACKAGE EVERYTHING FOR DOWNSTREAM ANALYSIS -- THREE TRIM LEVELS
% LEVEL 1 (sync_trim)     : raw hysteresis low-high-low boundaries, no buffer.
% LEVEL 2 (movement_trim) : Level 1 pulled in by movement_buffer_sec (default 4 s) each end.
% LEVEL 3 (trial_data)    : Level 2 with the first/last physical swing cycle
%                           removed via zero-crossing (this is the original
%                           "final" output -- kept under the same field
%                           names/filename as before for backward compatibility).
if ~exist('csv_out_dir', 'var')
    csv_out_dir = pwd;
end
% imu_fname_tag = matlab.lang.makeValidName(imu_file);

% Strip extension (.txt etc.) so the tag is e.g. KUKA_neckEMG_P06_Trial_Block1_01
% instead of ..._01_txt from makeValidName on the full filename.
[~, imu_base, ~] = fileparts(imu_file);
imu_fname_tag = matlab.lang.makeValidName(imu_base);

%% 9a. DIAGNOSTIC TRIM PLOTS (IMU + EMG, full recording with every trim
% boundary marked: Level-1 sync trim, Level-2 movement buffer (sync +/-
% movement_buffer_sec), and the Level-3 physical-cycle trim from Section
% 6/7). Purely diagnostic -- wrapped in try/catch so a plotting problem
% never takes down the batch run. Disable via P.make_diagnostic_plots = false.
if ~isfield(P, 'make_diagnostic_plots') || P.make_diagnostic_plots
    try
        %% --- IMU diagnostic plot ---
        % Sync channel: the raw signal the extractor thresholded on
        % (sync_extract_test_aug22.m), over the WHOLE recording -- same
        % 'Sync' column / last-column fallback that extractor uses.
        if ismember('Sync', imu_data.Properties.VariableNames)
            imu_sync_full = imu_data.Sync;
        else
            imu_sync_arr  = table2array(imu_data);
            imu_sync_full = imu_sync_arr(:, end);
        end

        % Main axis: rebuild the SAME direction-projected, rotation_lp_cutoff_hz
        % -lowpassed flange signal used for cycle detection in Section 5/6, but
        % over the full untrimmed imu_cal (not imu_cal_trimmed), so the sync,
        % buffer, and cycle boundaries can all be shown on one continuous
        % trace instead of three separately-scaled ones. Reuses imu_cal_lp1
        % (the full-recording lowpass computed once in Section 3, ahead of
        % any trimming -- see the note there) and coeff (for the Dia case)
        % exactly as computed in Section 5, instead of re-running filtfilt
        % here. No mean subtraction, same reasoning as Section 5: gyro DC
        % offset is already removed by MEMS calibration.
        g_flange_lp_fullrec = imu_cal_lp1(:, FLANGE_IMU_COLS(4:6));
        g_flange_xy_fullrec = g_flange_lp_fullrec(:, 1:2);
        if strcmp(direction, 'ML')
            flange_rot_fullrec = g_flange_xy_fullrec(:, 1);
        elseif strcmp(direction, 'AP')
            flange_rot_fullrec = g_flange_xy_fullrec(:, 2);
        elseif contains(direction, 'Dia')
            flange_rot_fullrec_proj = g_flange_xy_fullrec * coeff;
            flange_rot_fullrec      = flange_rot_fullrec_proj(:, 1);
        else
            flange_rot_fullrec = g_flange_xy_fullrec(:, 1);  % Section 5 would already have errored before reaching here
        end
        t_imu_fullrec = (0:length(flange_rot_fullrec)-1) / imu_fs;

        % idxStart_imu/idxEnd_imu (Section 2) and idxStart_imu_mv/idxEnd_imu_mv
        % (Section 2b) already index into the full recording. idxStart_imu_cycle/
        % idxEnd_imu_cycle (Section 6) are relative to the sync-trimmed array --
        % shift by idxStart_imu to land back on the full-recording index.
        abs_idxStart_imu_cycle = idxStart_imu + idxStart_imu_cycle - 1;
        abs_idxEnd_imu_cycle   = idxStart_imu + idxEnd_imu_cycle   - 1;

        imu_marks = struct( ...
            'sync_start',   idxStart_imu            / imu_fs, ...
            'sync_end',     idxEnd_imu              / imu_fs, ...
            'buffer_start', idxStart_imu_mv         / imu_fs, ...
            'buffer_end',   idxEnd_imu_mv           / imu_fs, ...
            'cycle_start',  abs_idxStart_imu_cycle  / imu_fs, ...
            'cycle_end',    abs_idxEnd_imu_cycle    / imu_fs);

        imu_plot_title = sprintf('%s -- Muscle %d, %s, %s (IMU)', imu_base, Muscle, direction, char(ParameterSet));
        imu_plot_path  = fullfile(csv_out_dir, sprintf('trimdiag_IMU_%s_Muscle%d_%s_%s.png', ...
            imu_fname_tag, Muscle, direction, matlab.lang.makeValidName(char(ParameterSet))));

        plot_trim_diagnostic(t_imu_fullrec, imu_sync_full, 'Sync channel', ...
            flange_rot_fullrec, sprintf('Flange %s axis (%.0f Hz LP)', direction, rotation_lp_cutoff_hz), ...
            imu_marks, imu_plot_title, imu_plot_path);

        %% --- EMG diagnostic plot ---
        % Sync channel: same raw channel the extractor thresholded on
        % (P.emg_sync_ch, e.g. 129), over the whole recording.
        if isfield(P, 'emg_sync_ch') && ~isempty(P.emg_sync_ch)
            emg_sync_full = double(stream.samples(P.emg_sync_ch, :));
        else
            emg_sync_full = nan(1, size(hp_emg_array, 2));  % keep the 2-panel layout even without a sync trace
            warning('P.emg_sync_ch not supplied -- EMG diagnostic plot will omit the sync trace.');
        end

        emg_main_full = hp_emg_array(1, :);   % channel 1 of the selected muscle, full recording, >50 Hz high-passed
        t_emg_fullrec = (0:size(hp_emg_array, 2)-1) / emg_fs;

        % idxStart_emg/idxEnd_emg (Section 2) and idxStart_emg_mv/idxEnd_emg_mv
        % (Section 2b) already index into the full recording. idxStart_emg_cycle/
        % idxEnd_emg_cycle (Section 7) are relative to the sync-trimmed array.
        abs_idxStart_emg_cycle = idxStart_emg + idxStart_emg_cycle - 1;
        abs_idxEnd_emg_cycle   = idxStart_emg + idxEnd_emg_cycle   - 1;

        emg_marks = struct( ...
            'sync_start',   idxStart_emg            / emg_fs, ...
            'sync_end',     idxEnd_emg              / emg_fs, ...
            'buffer_start', idxStart_emg_mv         / emg_fs, ...
            'buffer_end',   idxEnd_emg_mv           / emg_fs, ...
            'cycle_start',  abs_idxStart_emg_cycle  / emg_fs, ...
            'cycle_end',    abs_idxEnd_emg_cycle    / emg_fs);

        emg_plot_title = sprintf('%s -- Muscle %d, %s, %s (EMG ch. %d)', ...
            imu_base, Muscle, direction, char(ParameterSet), ch_range(1));
        emg_plot_path  = fullfile(csv_out_dir, sprintf('trimdiag_EMG_%s_Muscle%d_%s_%s.png', ...
            imu_fname_tag, Muscle, direction, matlab.lang.makeValidName(char(ParameterSet))));

        plot_trim_diagnostic(t_emg_fullrec, emg_sync_full, 'Sync channel', ...
            emg_main_full, sprintf('EMG ch. %d (>50 Hz HP)', ch_range(1)), ...
            emg_marks, emg_plot_title, emg_plot_path);

        fprintf('Saved trim diagnostic plots:\n  %s\n  %s\n', imu_plot_path, emg_plot_path);
    catch ME
        warning('Diagnostic trim plotting failed (non-fatal): %s', ME.message);
    end
end

% % --- LEVEL 1: sync trim (full Level-1-length arrays, no slicing needed)

% data_sync = package_level_dataset(emg_bp_all_sync, emg_rect_all_sync, ch_range, Muscle, f_low, f_high, emg_fs, ...
%     flange_rot, head_rot, torso_rot, flange_rot_broadband, head_rot_broadband, ...
%     imu_cal_trimmed(:, FLANGE_IMU_COLS), imu_cal_trimmed(:, HEAD_IMU_COLS), ...
%     imu_cal_trimmed(:, TORSO_IMU_COLS(1:6*double(num_imus_detected==3))), imu_fs, num_imus_detected);
% data_sync.meta = build_meta(direction, mu, T_cycle, num_cycles, idxStart_emg, idxEnd_emg, ...
%     idxStart_imu, idxEnd_imu, num_imus_detected, imu_file, rec_num, movement_buffer_sec, 'sync_trim');
% sync_mat_path = fullfile(csv_out_dir, sprintf('sync_trim_%s_Muscle%d.mat', imu_fname_tag, Muscle));
% data_sync.meta.mat_path = sync_mat_path;
% save(sync_mat_path, 'data_sync', '-v7.3');
% fprintf('Saved Level 1 (sync trim) data to:\n  %s\n', sync_mat_path);

% % --- LEVEL 2: movement trim (slice Level-1 arrays at the movement buffer) ---
% data_movement = package_level_dataset( ...
%     emg_bp_all_sync(:, mv_start_rel_emg:mv_end_rel_emg), ...
%     emg_rect_all_sync(:, mv_start_rel_emg:mv_end_rel_emg), ch_range, Muscle, f_low, f_high, emg_fs, ...
%     flange_rot(mv_start_rel_imu:mv_end_rel_imu), head_rot(mv_start_rel_imu:mv_end_rel_imu), ...
%     torso_rot(mv_start_rel_imu:mv_end_rel_imu), ...
%     flange_rot_broadband(mv_start_rel_imu:mv_end_rel_imu), head_rot_broadband(mv_start_rel_imu:mv_end_rel_imu), ...
%     imu_cal_trimmed(mv_start_rel_imu:mv_end_rel_imu, FLANGE_IMU_COLS), ...
%     imu_cal_trimmed(mv_start_rel_imu:mv_end_rel_imu, HEAD_IMU_COLS), ...
%     imu_cal_trimmed(mv_start_rel_imu:mv_end_rel_imu, TORSO_IMU_COLS(1:6*double(num_imus_detected==3))), ...
%     imu_fs, num_imus_detected);
% data_movement.meta = build_meta(direction, mu, T_cycle, num_cycles, idxStart_emg_mv, idxEnd_emg_mv, ...
%     idxStart_imu_mv, idxEnd_imu_mv, num_imus_detected, imu_file, rec_num, movement_buffer_sec, 'movement_trim');
% movement_mat_path = fullfile(csv_out_dir, sprintf('movement_trim_%s_Muscle%d.mat', imu_fname_tag, Muscle));
% data_movement.meta.mat_path = movement_mat_path;
% save(movement_mat_path, 'data_movement', '-v7.3');
% fprintf('Saved Level 2 (movement trim, %.1f s buffer) data to:\n  %s\n', movement_buffer_sec, movement_mat_path);

% % --- LEVEL 3: trial data (cycle-trim -- the original main output) ---
% data = struct();
% data.status = 'ok';
% 
% data.emg.signal        = EMG_all;                 % [n_emg_channels x N_final]
% data.emg.t              = t_emg_final;             % re-zeroed time base (s), matches IMU resampled signals
% data.emg.fs             = emg_fs;
% data.emg.channel_ids_global = ch_range;             % global channel numbers (1-128) matching row order
% data.emg.muscle         = Muscle;
% data.emg.bp_f_low       = f_low;
% data.emg.bp_f_high      = f_high;
% data.emg.signal_broadband_rect = EMG_broadband_rect_all; % [n_emg_channels x N_final], same length/time base as data.emg.signal; high-pass(>50Hz)+rectified ONLY (no mu-centered bandpass) -- use for time-frequency/spectral analyses, NOT for CRP
% 
% data.imu_flange.rotation_native   = flange_rot_cycle;      % cycle-trimmed, native imu_fs
% data.imu_flange.rotation_resamp   = Flange_rot_2000Hz;      % resampled onto data.emg.t / N_final
% data.imu_flange.six_axis_native   = imu_flange_6axis_cycle; % [N_imu_cycle x 6] = [Ax Ay Az Gx Gy Gz], calibrated, cycle-trimmed
% data.imu_flange.t_native          = t_imu_cycle;
% data.imu_flange.fs_native         = imu_fs;
% data.imu_flange.rotation_native_broadband = flange_rot_broadband_cycle; % same axis/PCA projection as rotation_native, but DC-removed only (no 1 Hz lowpass) -- use for time-frequency plots
% 
% data.imu_head.rotation_native   = head_rot_cycle;
% data.imu_head.rotation_resamp   = Head_rot_2000Hz;
% data.imu_head.six_axis_native   = imu_head_6axis_cycle;
% data.imu_head.t_native          = t_imu_cycle;
% data.imu_head.fs_native         = imu_fs;
% data.imu_head.rotation_native_broadband = head_rot_broadband_cycle; % same axis/PCA projection as rotation_native, but DC-removed only (no 1 Hz lowpass) -- use for time-frequency plots
% 
% if num_imus_detected == 3
%     data.imu_torso.rotation_native = torso_rot_cycle;
%     data.imu_torso.rotation_resamp = Torso_rot_2000Hz;
%     data.imu_torso.six_axis_native = imu_torso_6axis_cycle;
%     data.imu_torso.t_native        = t_imu_cycle;
%     data.imu_torso.fs_native       = imu_fs;
% end
% 
% data.meta = build_meta(direction, mu, T_cycle, num_cycles, idxStart_emg_cycle, idxEnd_emg_cycle, ...
%     idxStart_imu_cycle, idxEnd_imu_cycle, num_imus_detected, imu_file, rec_num, movement_buffer_sec, 'trial_data');
% data.meta.dur_sync_diff_pct   = dur_diff_pct;
% data.meta.dur_cycle_diff_pct  = dur_cycle_diff_pct;
% % Absolute (whole-recording) index versions, since idxStart/EndIdx_*_cycle
% % above are relative to the movement-trimmed array, not the raw recording.
% data.meta.idxStart_imu_cycle_abs = idxStart_imu_mv + idxStart_imu_cycle - mv_start_rel_imu;
% data.meta.idxEnd_imu_cycle_abs   = idxStart_imu_mv + idxEnd_imu_cycle   - mv_start_rel_imu;
% data.meta.idxStart_emg_cycle_abs = idxStart_emg_mv + idxStart_emg_cycle - mv_start_rel_emg;
% data.meta.idxEnd_emg_cycle_abs   = idxStart_emg_mv + idxEnd_emg_cycle   - mv_start_rel_emg;
% % data.meta.sync_trim_mat_path     = sync_mat_path;
% % data.meta.movement_trim_mat_path = movement_mat_path;
% 
% data_mat_path = fullfile(csv_out_dir, sprintf('analysis_data_%s_Muscle%d.mat', imu_fname_tag, Muscle));
% data.meta.mat_path = data_mat_path;   % so callers/batch driver know where this pair's file landed
% save(data_mat_path, 'data', '-v7.3');
% fprintf('Saved Level 3 (trial data, cycle-trimmed) data (%d EMG channels, flange + head IMU) to:\n  %s\n\n', ...
%     n_emg_channels, data_mat_path);

% --- LEVEL 3: trial data (cycle-trim -- the original main output) ---

% --- Sync-trimmed-only EMG (former Level 1 output, folded in here) ---
% Sync-trim length only (idxStart_emg:idxEnd_emg) -- NOT cropped to the
% physical-swing cycle boundaries like data.emg.signal / signal_broadband_rect
% above, so these are a different length (N_sync, not N_final) and carry
% their own time base, data.emg.t_synctrim.


data.emg.signal        = EMG_all;                 % [n_emg_channels x N_final]
data.emg.t              = t_emg_final;             % re-zeroed time base (s), matches IMU resampled signals
data.emg.fs             = emg_fs;
data.emg.channel_ids_global = ch_range;             % global channel numbers (1-128) matching row order
data.emg.muscle         = Muscle;
data.emg.bp_f_low       = f_low;
data.emg.bp_f_high      = f_high;
data.emg.signal_broadband_rect = EMG_broadband_rect_all; % [n_emg_channels x N_final], same length/time base as data.emg.signal; high-pass(>50Hz)+rectified ONLY (no mu-centered bandpass) -- use for time-frequency/spectral analyses, NOT for CRP
data.emg.signal_rms            = EMG_rms_all;            % [n_emg_channels x N_final], same length/time base as data.emg.signal; high-pass(>50Hz) + centered RMS (rms_win_sec window, default 400ms) ONLY (no mu-centered bandpass, no rectify) -- cycle-trimmed, same idxStart/idxEnd_emg_cycle boundaries as signal_broadband_rect
data.emg.signal_highpass_cycletrim = EMG_highpass_cycle_all; % [n_emg_channels x N_final], same length/time base as data.emg.signal; high-pass(>50Hz) ONLY -- not rectified, not RMS'd, not bandpassed; cycle-trimmed to the same boundaries as signal_broadband_rect/signal_rms
data.emg.rms_win_sec            = emg_rms_window_sec;    % window length (s) used for signal_rms / signal_synctrim_rms (centered/non-causal)

% --- Per-channel d-prime (Section 8d) ---
% Gaussian-band-filtered (centered at mu), forward-lag-aligned (lag > 0
% only, lag = 0 excluded), velocity-binned d-prime sensitivity curve for
% every EMG channel. See Section 8d above for the full method.
data.dprime.values             = dprime_all;              % [n_emg_channels x num_bins], d' per channel per velocity bin
data.dprime.bin_centers        = dprime_bin_centers;       % [1 x num_bins], deg/s
data.dprime.off_bin_idx        = dprime_off_bin_idx;        % index into bin_centers/values used as the near-0 "off" baseline
data.dprime.lag_sec            = dprime_lags_sec;           % [n_emg_channels x 1], optimal forward (>0) lag per channel, s
data.dprime.channel_ids_global = ch_range;                  % global channel numbers (1-128), matches data.emg.channel_ids_global row order
data.dprime.gaussian_fwhm_hz   = dprime_fwhm;
data.dprime.gaussian_center_hz = mu;
data.dprime.bin_width_deg_s    = dprime_bin_width;
data.dprime.xcorr_max_lag_sec  = dprime_max_lag_sec;
data.dprime.min_samples_per_bin = dprime_min_samples_per_bin;

% --- Sync-trimmed-only EMG (former Level 1 output, folded in here) ---
% Sync-trim length only (idxStart_emg:idxEnd_emg) -- NOT cropped to the
% physical-swing cycle boundaries like data.emg.signal / signal_broadband_rect
% above, so these are a different length (N_sync, not N_final) and carry
% their own time base, data.emg.t_synctrim.
data.emg.signal_synctrim_highpass = hp_emg_trimmed;     % [n_emg_channels x N_sync], sync-trimmed, high-pass ONLY -- not rectified, not bandpassed
data.emg.signal_synctrim_rect     = emg_rect_all_sync;   % [n_emg_channels x N_sync], sync-trimmed, high-pass + rectified
data.emg.signal_synctrim_rms      = emg_rms_all_sync;    % [n_emg_channels x N_sync], sync-trimmed, high-pass + centered RMS (rms_win_sec window) -- NOT cycle-trimmed
data.emg.t_synctrim               = t_emg_trimmed;       % time base (s) for the three arrays above


% Flange IMU
data.imu_flange.rotation_native   = flange_rot_cycle;      
data.imu_flange.rotation_resamp   = Flange_rot_2000Hz;      
data.imu_flange.six_axis_native   = imu_flange_6axis_cycle; % [N_imu_cycle x 6] @ native fs
data.imu_flange.six_axis_resamp   = Flange_6axis_2000Hz;   % [N_final x 6] @ emg_fs, 5 Hz lowpass (filtfilt) applied at native imu_fs before resampling
data.imu_flange.t_native          = t_imu_cycle;
data.imu_flange.fs_native         = imu_fs;
data.imu_flange.rotation_native_broadband = flange_rot_broadband_cycle;
% Sync-trimmed (NOT cycle-trimmed), resampled onto the sync-trimmed EMG
% time base (N_sync_only samples @ emg_fs) -- requirement: sync-trimmed
% IMU, resampled. Same idea as rotation_resamp/six_axis_resamp above, but
% for the full Level-1 sync-trim length instead of the cycle-trimmed one.
data.imu_flange.rotation_synctrim_resamp  = Flange_rot_sync_only_resamp;  % [N_sync_only x 1]
data.imu_flange.six_axis_synctrim_resamp  = Flange_6axis_synctrim_resamp; % [N_sync_only x 6], 5 Hz lowpass applied at native imu_fs before resampling
data.imu_flange.t_synctrim_resamp         = t_sync_only;                  % time base (s) for the two arrays above; numerically identical to data.emg.t_synctrim (same N_sync length @ emg_fs)

% Head IMU
data.imu_head.rotation_native   = head_rot_cycle;
data.imu_head.rotation_resamp   = Head_rot_2000Hz;
data.imu_head.six_axis_native   = imu_head_6axis_cycle;  % [N_imu_cycle x 6] @ native fs
data.imu_head.six_axis_resamp   = Head_6axis_2000Hz;     % [N_final x 6] @ emg_fs, 5 Hz lowpass (filtfilt) applied at native imu_fs before resampling
data.imu_head.t_native          = t_imu_cycle;
data.imu_head.fs_native         = imu_fs;
data.imu_head.rotation_native_broadband = head_rot_broadband_cycle;
data.imu_head.rotation_synctrim_resamp  = Head_rot_sync_only_resamp;   % [N_sync_only x 1]
data.imu_head.six_axis_synctrim_resamp  = Head_6axis_synctrim_resamp;  % [N_sync_only x 6], 5 Hz lowpass applied at native imu_fs before resampling
data.imu_head.t_synctrim_resamp         = t_sync_only;

% Torso IMU (if present)
if num_imus_detected == 3
    data.imu_torso.rotation_native = torso_rot_cycle;
    data.imu_torso.rotation_resamp = Torso_rot_2000Hz;
    data.imu_torso.six_axis_native = imu_torso_6axis_cycle; % [N_imu_cycle x 6] @ native fs
    data.imu_torso.six_axis_resamp = Torso_6axis_2000Hz;   % [N_final x 6] @ emg_fs, 5 Hz lowpass (filtfilt) applied at native imu_fs before resampling
    data.imu_torso.t_native        = t_imu_cycle;
    data.imu_torso.fs_native       = imu_fs;
    data.imu_torso.rotation_synctrim_resamp = Torso_rot_sync_only_resamp;  % [N_sync_only x 1]
    data.imu_torso.six_axis_synctrim_resamp = Torso_6axis_synctrim_resamp; % [N_sync_only x 6]
    data.imu_torso.t_synctrim_resamp        = t_sync_only;
end

data.meta = build_meta(direction, mu, T_cycle, num_cycles, idxStart_emg_cycle, idxEnd_emg_cycle, ...
    idxStart_imu_cycle, idxEnd_imu_cycle, num_imus_detected, imu_file, rec_num, movement_buffer_sec, 'trial_data');
data.meta.dur_sync_diff_pct   = dur_diff_pct;
data.meta.dur_cycle_diff_pct  = dur_cycle_diff_pct;
data.meta.imu_6axis_lp_cutoff_hz = imu_6axis_lp_cutoff_hz;  % lowpass applied to *.six_axis_resamp before resampling to emg_fs
data.meta.rotation_lp_cutoff_hz  = rotation_lp_cutoff_hz;   % lowpass applied to the direction-projected rotation signal (cycle detection / *.rotation_resamp)
% Absolute (whole-recording) index versions, since idxStart/EndIdx_*_cycle
% above are relative to the movement-trimmed array, not the raw recording.
data.meta.idxStart_imu_cycle_abs = idxStart_imu_mv + idxStart_imu_cycle - mv_start_rel_imu;
data.meta.idxEnd_imu_cycle_abs   = idxStart_imu_mv + idxEnd_imu_cycle   - mv_start_rel_imu;
data.meta.idxStart_emg_cycle_abs = idxStart_emg_mv + idxStart_emg_cycle - mv_start_rel_emg;
data.meta.idxEnd_emg_cycle_abs   = idxStart_emg_mv + idxEnd_emg_cycle   - mv_start_rel_emg;
% data.meta.sync_trim_mat_path     = sync_mat_path;
% data.meta.movement_trim_mat_path = movement_mat_path;

% Filename: trimmed_<imu_base>_Muscle<N>_<direction>_<ParameterSet>_<condition>.mat
% e.g. trimmed_KUKA_neckEMG_P06_Trial_Block1_01_Muscle1_AP_1800_0504_EO.mat
data_mat_path = fullfile(csv_out_dir, sprintf('trimmed_%s_Muscle%d_%s_%s_%s.mat', ...
    imu_fname_tag, Muscle, direction, ParameterSet, condition));
data.meta.mat_path = data_mat_path;   % so callers/batch driver know where this pair's file landed
save(data_mat_path, 'data', '-v7.3');
fprintf('Saved Level 3 (trial data, cycle-trimmed) data (%d EMG channels, flange + head IMU) to:\n  %s\n\n', ...
    n_emg_channels, data_mat_path);

end % of emg_imu_analysis_data_prep_v2

%% ========================================================================
%  LEVEL-PACKAGING HELPERS (shared by all three trim levels in Section 9)
% ========================================================================

function level = package_level_dataset(emg_bp_slice, emg_rect_slice, ch_range, Muscle, f_low, f_high, emg_fs, ...
    flange_rot_slice, head_rot_slice, torso_rot_slice, flange_rot_bb_slice, head_rot_bb_slice, ...
    imu_flange_6ax_slice, imu_head_6ax_slice, imu_torso_6ax_slice, imu_fs, num_imus_detected)
% Builds one level's EMG + IMU (flange/head/torso) output struct, given
% already-sliced arrays for that level. Resamples the IMU rotation
% signals onto that level's own EMG-length time base, same as Section 8
% does for the final (cycle-trimmed) level.

    N_emg = size(emg_bp_slice, 2);
    N_imu = length(flange_rot_slice);
    t_imu_native = (0:N_imu-1) / imu_fs;

    level = struct();
    level.emg.signal                = emg_bp_slice;
    level.emg.signal_broadband_rect = emg_rect_slice;
    level.emg.t                     = (0:N_emg-1) / emg_fs;
    level.emg.fs                    = emg_fs;
    level.emg.channel_ids_global    = ch_range;
    level.emg.muscle                = Muscle;
    level.emg.bp_f_low              = f_low;
    level.emg.bp_f_high             = f_high;

    level.imu_flange.rotation_native            = flange_rot_slice;
    level.imu_flange.rotation_native_broadband  = flange_rot_bb_slice;
    level.imu_flange.rotation_resamp            = resample_to_length(flange_rot_slice, N_emg);
    level.imu_flange.six_axis_native            = imu_flange_6ax_slice;
    level.imu_flange.t_native                   = t_imu_native;
    level.imu_flange.fs_native                  = imu_fs;

    level.imu_head.rotation_native           = head_rot_slice;
    level.imu_head.rotation_native_broadband = head_rot_bb_slice;
    level.imu_head.rotation_resamp           = resample_to_length(head_rot_slice, N_emg);
    level.imu_head.six_axis_native           = imu_head_6ax_slice;
    level.imu_head.t_native                  = t_imu_native;
    level.imu_head.fs_native                 = imu_fs;

    if num_imus_detected == 3
        level.imu_torso.rotation_native = torso_rot_slice;
        level.imu_torso.rotation_resamp = resample_to_length(torso_rot_slice, N_emg);
        level.imu_torso.six_axis_native = imu_torso_6ax_slice;
        level.imu_torso.t_native        = t_imu_native;
        level.imu_torso.fs_native       = imu_fs;
    end
end

function resamp = resample_to_length(sig, N_target)
    sig = sig(:);
    resamp = interp1(linspace(0,1,length(sig)), sig, linspace(0,1,N_target), 'linear')';
end

function meta = build_meta(direction, mu, T_cycle, num_cycles, idxStart_emg, idxEnd_emg, ...
    idxStart_imu, idxEnd_imu, num_imus_detected, imu_file, rec_num, movement_buffer_sec, level_name)
    meta = struct();
    meta.level                 = level_name;
    meta.direction              = direction;
    meta.mu_hz                  = mu;
    meta.T_cycle_sec            = T_cycle;
    meta.num_cycles_detected    = num_cycles;
    meta.idxStart_emg           = idxStart_emg;   meta.idxEnd_emg = idxEnd_emg;
    meta.idxStart_imu           = idxStart_imu;   meta.idxEnd_imu = idxEnd_imu;
    meta.num_imus_detected      = num_imus_detected;
    meta.imu_file                = imu_file;
    meta.rec_num                 = rec_num;
    meta.movement_buffer_sec     = movement_buffer_sec;
end

%% ========================================================================
%  SUPPORTING HELPER FUNCTIONS (same framework as emg_imu_sync_pipeline_FIXED.m)
% ========================================================================

% extract_sync_trial_hysteresis / extract_all_segments / identify_trials_idx
% (the duplicated hysteresis low-high-low sync detector) have been
% removed -- sync-trim indices now come from sync_extract_test_aug22.m's
% EMG_IMU_Matched_Trials.csv via P (see Section 2 above), so this file no
% longer re-derives them from raw signals.


function plot_trim_diagnostic(t_vec, sig_top, top_label, sig_bottom, bottom_label, marks, title_str, out_png)
% Two-panel (sync channel on top, main axis on bottom) diagnostic figure,
% linked x-axis, with vertical lines marking the sync trim, movement
% buffer trim, and physical-cycle trim boundaries (all in seconds, via
% `marks.sync_start/sync_end/buffer_start/buffer_end/cycle_start/cycle_end`).
% Any mark left NaN is simply skipped. Saved to out_png and closed
% (Visible off throughout) so this is safe to call in a large batch loop.
    max_pts = 20000;
    n = numel(t_vec);
    if n > max_pts
        step     = ceil(n / max_pts);
        t_plot   = t_vec(1:step:end);
        top_plot = sig_top(1:step:end);
        bot_plot = sig_bottom(1:step:end);
    else
        t_plot = t_vec; top_plot = sig_top; bot_plot = sig_bottom;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 700]);

    ax1 = subplot(2, 1, 1);
    plot(ax1, t_plot, top_plot, 'Color', [0.4 0.4 0.4]);
    ylabel(ax1, top_label);
    title(ax1, title_str, 'Interpreter', 'none');
    add_trim_lines(ax1, marks);

    ax2 = subplot(2, 1, 2);
    plot(ax2, t_plot, bot_plot, 'Color', [0 0.35 0.65]);
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, bottom_label);
    add_trim_lines(ax2, marks);

    % Dummy lines purely to build a legend (xline handles are hidden from
    % the legend via HandleVisibility off inside add_trim_lines).
    hold(ax2, 'on');
    h_sync = plot(ax2, NaN, NaN, '-',  'Color', [0.85 0.10 0.10], 'LineWidth', 1.5);
    h_buf  = plot(ax2, NaN, NaN, '--', 'Color', [0.90 0.60 0.10], 'LineWidth', 1.3);
    h_cyc  = plot(ax2, NaN, NaN, '-.', 'Color', [0.10 0.60 0.20], 'LineWidth', 1.3);
    hold(ax2, 'off');
    legend(ax2, [h_sync, h_buf, h_cyc], {'Sync trim', 'Buffer trim (\pm movement buffer)', 'Cycle trim'}, ...
        'Location', 'best');

    linkaxes([ax1, ax2], 'x');

    print(fig, out_png, '-dpng', '-r100');
    close(fig);
end

function add_trim_lines(ax, m)
    if isfield(m, 'sync_start') && ~isnan(m.sync_start)
        xline(ax, m.sync_start, '-',  'Color', [0.85 0.10 0.10], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        xline(ax, m.sync_end,   '-',  'Color', [0.85 0.10 0.10], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    if isfield(m, 'buffer_start') && ~isnan(m.buffer_start)
        xline(ax, m.buffer_start, '--', 'Color', [0.90 0.60 0.10], 'LineWidth', 1.3, 'HandleVisibility', 'off');
        xline(ax, m.buffer_end,   '--', 'Color', [0.90 0.60 0.10], 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end
    if isfield(m, 'cycle_start') && ~isnan(m.cycle_start)
        xline(ax, m.cycle_start, '-.', 'Color', [0.10 0.60 0.20], 'LineWidth', 1.3, 'HandleVisibility', 'off');
        xline(ax, m.cycle_end,   '-.', 'Color', [0.10 0.60 0.20], 'LineWidth', 1.3, 'HandleVisibility', 'off');
    end
end

function idx_cross = find_debounced_zero_cross(sig, win_start, win_end, debounce_len, mode)
% Finds stable zero-crossings (either direction) within a search window.
% mode: 'first' -> earliest crossing; 'last' -> latest crossing;
%       'all'   -> full sorted list of absolute indices (into sig).
    sub_sig = sig(win_start:win_end);

    rising_rel  = find(sub_sig(1:end-1) <= 0 & sub_sig(2:end) > 0);
    falling_rel = find(sub_sig(1:end-1) >= 0 & sub_sig(2:end) < 0);

    valid_crossings = [];

    for k = 1:length(rising_rel)
        rel_i = rising_rel(k);
        check_end = min(length(sub_sig), rel_i + debounce_len);
        if all(sub_sig(rel_i + 1 : check_end) > 0)
            valid_crossings(end+1) = win_start + rel_i; %#ok<AGROW>
        end
    end

    for k = 1:length(falling_rel)
        rel_i = falling_rel(k);
        check_end = min(length(sub_sig), rel_i + debounce_len);
        if all(sub_sig(rel_i + 1 : check_end) < 0)
            valid_crossings(end+1) = win_start + rel_i; %#ok<AGROW>
        end
    end

    valid_crossings = sort(valid_crossings);

    if isempty(valid_crossings)
        idx_cross = [];
    elseif strcmp(mode, 'first')
        idx_cross = valid_crossings(1);
    elseif strcmp(mode, 'last')
        idx_cross = valid_crossings(end);
    else % 'all'
        idx_cross = valid_crossings;
    end
end

function freq_hz = parse_freq_estimate_from_paramset(paramset_str)
% ParameterSet strings look like "1320_0300": the leading number, divided
% by 1000, is the swing-frequency estimate in Hz (e.g. "1320" -> 1.320 Hz).
% Returns NaN if the string can't be parsed so callers can fall back to
% the estimateNumCycles()-based guess instead.
    freq_hz = NaN;
    if isempty(paramset_str)
        return;
    end
    tok = regexp(char(paramset_str), '^(\d+)', 'tokens', 'once');
    if isempty(tok)
        return;
    end
    freq_hz = str2double(tok{1}) / 1000;
    if ~isfinite(freq_hz) || freq_hz <= 0
        freq_hz = NaN;
    end
end

function calibR = createCalibration(c1, c2)
    v1 = c1(:) / norm(c1);
    v2 = c2(:) / norm(c2);

    z_axis = v1;
    y_axis = cross(z_axis, v2);
    if norm(y_axis) < 1e-6
        y_axis = [0; 1; 0];
    else
        y_axis = y_axis / norm(y_axis);
    end
    x_axis = cross(y_axis, z_axis);

    calibR = [x_axis, y_axis, z_axis]';
end

% function calibrated_data = apply_mems_calibration(raw_data, calib_files)
%     b_gyro = load(calib_files.b_gyro);
%     b_acc  = load(calib_files.b_acc);
%     s_acc  = load(calib_files.s_acc);
%     t_acc  = load(calib_files.t_acc);
% 
%     b_gyro = b_gyro(:);
%     b_acc  = b_acc(:);
% 
%     A_raw = raw_data(:, 1:3)';
%     G_raw = raw_data(:, 4:6)';
% 
%     A_unbiased = A_raw - b_acc;
%     A_cal      = t_acc * (s_acc * A_unbiased);
%     G_cal      = G_raw - b_gyro;
% 
% 
%     calibrated_data = [A_cal', G_cal'];
% end

function calibrated_data = apply_mems_calibration(raw_data, calib_files)
    b_gyro = load(calib_files.b_gyro);
    b_acc  = load(calib_files.b_acc);
    s_acc  = load(calib_files.s_acc);
    t_acc  = load(calib_files.t_acc);

    b_gyro = b_gyro(:)';
    b_acc  = b_acc(:)';

    A_raw = raw_data(:, 1:3);
    G_raw = raw_data(:, 4:6);

    A_unbiased = A_raw - b_acc;
    A_cal      = t_acc * (s_acc * A_unbiased);
    G_cal      = G_raw - b_gyro;


    calibrated_data = [A_cal, G_cal];
end

function cal_imu = calibrateIMU(imu_raw, calibR)
    cal_imu = imu_raw;
    if size(imu_raw, 2) >= 6
        cal_imu(:, 1:3) = (calibR * imu_raw(:, 1:3)')';
        cal_imu(:, 4:6) = (calibR * imu_raw(:, 4:6)')';
    end
end

function numCycles = estimateNumCycles(signal)
    absSignal = abs(signal - mean(signal));

    maxVal = max(absSignal);
    [peaks, ~] = findpeaks(absSignal, 'MinPeakHeight', maxVal * 0.50);

    if isempty(peaks)
        peaks = maxVal;
    end

    sortedPeaks = sort(peaks, 'descend');
    numPeaksToAverage = min(5, length(sortedPeaks));
    avgMaxVal = mean(sortedPeaks(1:numPeaksToAverage));

    swingThreshold = avgMaxVal * 0.75;

    aboveThreshold = absSignal >= swingThreshold;
    numSwings = sum(diff(aboveThreshold) == 1);

    numCycles = numSwings / 2;
end

function [imu1_6ax, imu2_6ax, imu4_6ax] = extract_calib_arrays_6axis(tbl)
% Like extract_calib_arrays() below, but returns the full
% [Ax Ay Az Gx Gy Gz] block per IMU instead of just 3 accel columns, so
% apply_mems_calibration has what it needs. Uses the same positional
% layout as extract_calib_arrays (accel-first-of-6-per-IMU, anchored at
% col 1 for headerless 12-26 col files, anchored at col 10 for >=27 col
% files). This anchor is also used as the FALLBACK when named columns
% (e.g. one of the six per-IMU labels) can't be matched, so it must be
% chosen based on this table's actual width, not assumed to be the
% >=27-col layout -- a named-but-narrower file (e.g. a 2-IMU AllOri file)
% would otherwise index past the end of the table.
    vars = tbl.Properties.VariableNames;
    has_named_headers = any(contains(vars, 'Ax'));
    num_cols = width(tbl);

    if num_cols >= 27
        fb1 = 10:15; fb2 = 16:21; fb4 = 22:27;
    elseif num_cols >= 18
        fb1 = 1:6;   fb2 = 7:12;  fb4 = 13:18;
    elseif num_cols >= 12
        fb1 = 1:6;   fb2 = 7:12;  fb4 = 1:6;   % no room for a 3rd IMU -- mirror imu4=imu1 fallback
    else
        fb1 = []; fb2 = []; fb4 = [];  % too narrow for any fallback; only a full named match can work
    end

    if has_named_headers
        imu1_6ax = get_named_6axis_block(tbl, vars, {'Ax1','Ay1','Az1','Gx1','Gy1','Gz1'}, ...
                                                       {'Ax_1','Ay_1','Az_1','Gx_1','Gy_1','Gz_1'}, fb1, 'IMU1');
        imu2_6ax = get_named_6axis_block(tbl, vars, {'Ax2','Ay2','Az2','Gx2','Gy2','Gz2'}, ...
                                                       {'Ax_2','Ay_2','Az_2','Gx_2','Gy_2','Gz_2'}, fb2, 'IMU2');
        % Torso may be labeled '3' or '4' depending on file -- try both.
        if any(ismember({'Ax4','Ax_4'}, vars))
            imu4_6ax = get_named_6axis_block(tbl, vars, {'Ax4','Ay4','Az4','Gx4','Gy4','Gz4'}, ...
                                                           {'Ax_4','Ay_4','Az_4','Gx_4','Gy_4','Gz_4'}, fb4, 'IMU4');
        elseif any(ismember({'Ax3','Ax_3'}, vars))
            imu4_6ax = get_named_6axis_block(tbl, vars, {'Ax3','Ay3','Az3','Gx3','Gy3','Gz3'}, ...
                                                           {'Ax_3','Ay_3','Az_3','Gx_3','Gy_3','Gz_3'}, fb4, 'IMU4');
        else
            imu4_6ax = imu1_6ax;  % no 3rd IMU present -- match original acc4=acc1 fallback
        end
    else
        raw = table2array(tbl);
        if num_cols >= 12
            imu1_6ax = raw(:, fb1);
            imu2_6ax = raw(:, fb2);
            imu4_6ax = raw(:, fb4);
        else
            error(['extract_calib_arrays_6axis: only %d columns found -- not enough to ' ...
                   'extract a second IMU''s full 6-axis block (need >=12).'], num_cols);
        end
    end
end

function block6 = get_named_6axis_block(tbl, vars, names_no_underscore, names_underscore, fallback_cols, label)
    if all(ismember(names_no_underscore, vars))
        block6 = table2array(tbl(:, names_no_underscore));
    elseif all(ismember(names_underscore, vars))
        block6 = table2array(tbl(:, names_underscore));
    elseif ~isempty(fallback_cols) && max(fallback_cols) <= width(tbl)
        raw = table2array(tbl);
        block6 = raw(:, fallback_cols);
    else
        error(['get_named_6axis_block: could not find all 6 named columns for %s, and this ' ...
               'table (%d cols) is too narrow for the positional fallback (needs col %d). ' ...
               'Check this file''s header naming against diagnose_calib_columns.m.'], ...
               label, width(tbl), max([fallback_cols, 0]));
    end
end

function pose_accel_mean = calibrated_pose_accel_mean(imu_6ax_block, calib_files)
% Applies the SAME gyro-offset removal + offset/scale/misalignment
% correction the trial data gets (Section 3's IMU1_AandG/IMU2_AandG
% treatment), to one static-pose recording (allori1/allori2), then
% returns the mean of the resulting (corrected) accelerometer columns
% for use in createCalibration().
    % imu_6ax_block(:, 4:6) = imu_6ax_block(:, 4:6) - mean(imu_6ax_block(:, 4:6), 1, 'omitnan');
    calibrated = apply_mems_calibration(imu_6ax_block, calib_files);
    pose_accel_mean = mean(calibrated(:, 1:3), 1, 'omitnan');
end

function [acc1, acc2, acc4] = extract_calib_arrays(tbl)
% NOTE: superseded by extract_calib_arrays_6axis (above) for building
% calibR_1/calibR_2/calibR_4 in Section 3, since that fix needs the full
% 6-axis block (not just accel) to run apply_mems_calibration on each
% pose. Left here unused in case anything else still calls it.
    vars = tbl.Properties.VariableNames;
    has_named_headers = any(contains(vars, 'Ax'));

    if has_named_headers
        if ismember('Ax1', vars) && ismember('Ay1', vars) && ismember('Az1', vars)
            acc1 = [tbl.Ax1, tbl.Ay1, tbl.Az1];
        elseif ismember('Ax_1', vars) && ismember('Ay_1', vars) && ismember('Az_1', vars)
            acc1 = [tbl.Ax_1, tbl.Ay_1, tbl.Az_1];
        else
            raw = table2array(tbl);
            acc1 = raw(:, 10:12);
        end

        if ismember('Ax2', vars) && ismember('Ay2', vars) && ismember('Az2', vars)
            acc2 = [tbl.Ax2, tbl.Ay2, tbl.Az2];
        elseif ismember('Ax_2', vars) && ismember('Ay_2', vars) && ismember('Az_2', vars)
            acc2 = [tbl.Ax_2, tbl.Ay_2, tbl.Az_2];
        else
            raw = table2array(tbl);
            acc2 = raw(:, 16:18);
        end

        if ismember('Ax4', vars) && ismember('Ay4', vars) && ismember('Az4', vars)
            acc4 = [tbl.Ax4, tbl.Ay4, tbl.Az4];
        elseif ismember('Ax_4', vars) && ismember('Ay_4', vars) && ismember('Az_4', vars)
            acc4 = [tbl.Ax_4, tbl.Ay_4, tbl.Az_4];
        elseif ismember('Ax3', vars) && ismember('Ay3', vars) && ismember('Az3', vars)
            acc4 = [tbl.Ax3, tbl.Ay3, tbl.Az3];
        elseif ismember('Ax_3', vars) && ismember('Ay_3', vars) && ismember('Az_3', vars)
            acc4 = [tbl.Ax_3, tbl.Ay_3, tbl.Az_3];
        else
            raw = table2array(tbl);
            if size(raw, 2) >= 24
                acc4 = raw(:, 22:24);
            else
                acc4 = acc1;
            end
        end

    else
        raw = table2array(tbl);
        num_cols = size(raw, 2);

        if num_cols >= 27
            acc1 = raw(:, 10:12);
            acc2 = raw(:, 16:18);
            acc4 = raw(:, 22:24);
        elseif num_cols >= 9
            acc1 = raw(:, 1:3);
            acc2 = raw(:, 4:6);
            acc4 = raw(:, 7:9);
        elseif num_cols >= 6
            acc1 = raw(:, 1:3);
            acc2 = raw(:, 4:6);
            acc4 = acc1;
        else
            error('Headerless calibration table has insufficient columns (%d found).', num_cols);
        end
    end
end

function log_skipped_sync(out_dir, fname, emg_ok, imu_ok)
    skipped_csv_path = fullfile(out_dir, 'Skipped_Sync_Trials.csv');
    reason_str = sprintf('EMG_OK=%d, IMU_OK=%d', emg_ok, imu_ok);

    tbl = table(string(fname), string(reason_str), 'VariableNames', {'FileName', 'Reason'});
    if exist(skipped_csv_path, 'file')
        writetable(tbl, skipped_csv_path, 'WriteMode', 'append');
    else
        writetable(tbl, skipped_csv_path);
    end
end

% get_bit_volts removed -- it only ever converted the EMG sync channel's
% raw samples to mV for the (now-removed) in-function sync detector.
% Nothing else in this file reads bit_volts.