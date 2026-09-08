%% ========================================================================
%  EMG / IMU BATCH DATA PREP DRIVER
% ========================================================================
%  Bridges sync_pulse_extractor_v5_15th_pathfix_qwen.m and
%  _analysis_dataprep_aug22_SyncLowPassInterp_plot.m (now a function -- see that file's
%  header). Reads the sync extractor's EMG_IMU_Match_Key.csv for one
%  participant and runs the data-prep pipeline once per matched pair,
%  automatically, instead of you setting rec_path/rec_num/imu_file by
%  hand for a single trial each time.
%
%  WORKFLOW:
%    1. Run the sync pulse extractor for this participant first. It
%       writes <csv_out_dir>/EMG_IMU_Match_Key.csv, one row per matched
%       EMG<->IMU trial pair (EMG_SessionIdx, EMG_LocalRec, IMU_FileName,
%       EMG_FirstHighToLowIdx, EMG_FinalLowToHighIdx, IMU_FirstHighToLowIdx,
%       IMU_FinalLowToHighIdx, ...) -- the extractor now folds those four
%       sync sample indices into the Match Key directly, so this is the
%       only file the driver needs to read for pairing + sync indices.
%    2. Run ClaudeKeySetting.m for this same participant. It writes
%       <key_setting_out_dir>/IMU_Key_Table.csv, one row per IMU file
%       (FileName, ParameterSet, Condition, GroupID, GroupLabel), where
%       GroupLabel is one of 'ML' / 'AP' / 'DiaR' / 'DiaL' -- exactly the
%       values _analysis_dataprep_aug22_SyncLowPassInterp_plot's `direction` input expects
%       (it does strcmp for 'ML'/'AP' and contains(...,'Dia') for
%       'DiaR'/'DiaL'). This script joins that table onto the Match Key
%       by IMU_FileName so direction no longer has to be set by hand --
%       it's read per-trial straight from the key table.
%    3. Set CONFIG below to match that same participant (rec_paths /
%       imu_folder / calib1_file / calib2_file MUST be identical to what
%       you used in the sync extractor, since EMG_SessionIdx indexes into
%       rec_paths and IMU_FileName must resolve inside imu_folder).
%    4. Set the processing-specific fields (MEMS calibration file structs)
%       that _analysis_dataprep_aug22_SyncLowPassInterp_plot needs but the sync extractor
%       doesn't. Muscle is no longer a single value here -- see
%       muscle_list below.
%    5. Run this script. For every row of the Key CSV, it looks up that
%       row's direction from the key table, then calls
%       _analysis_dataprep_aug22_SyncLowPassInterp_plot ONCE PER ENTRY OF muscle_list
%       (default [1 2], i.e. both electrode banks / all 128 channels
%       across the two calls), writing one .mat file per
%       (pair, Muscle) combination (analysis_data_<imu_file>_Muscle<N>.mat)
%       plus a single Batch_Data_Prep_Log.csv summarizing every run.
%       One (pair, Muscle) combination failing (error, sync-validation
%       skip, or missing direction) does NOT stop the batch -- it's
%       logged and the loop moves on.
% ========================================================================

clear; clc;

%% ------------------------- CONFIG --------------------------------------
participant_num = 10;  % must match the participant_num used in the sync extractor


% --- Must match the sync_pulse_extractor CONFIG for this participant ---
if participant_num == 1
    rec_paths   = {'/media/veracrypt1/P01/2026-07-24_14-02-02/'};
    imu_folder  = '/media/veracrypt1/P01';
    calib1_file = 'KUKA_neckEMG_P01_Trial_calib1_redo2.txt';
    calib2_file = 'KUKA_neckEMG_P01_Trial_calib2.txt';
elseif participant_num == 2
    rec_paths = {
        '/media/veracrypt1/P02/2026-07-27_10-15-17/', ...
        '/media/veracrypt1/P02/2026-07-27_12-03-56/'
    };
    imu_folder  = '/media/veracrypt1/P02';
    calib1_file = 'KUKA_neckEMG_P02_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P02_AllOri2.txt';
elseif participant_num == 3
    rec_paths   = {'/media/veracrypt1/P03/2026-07-27_14-37-14/'};
    imu_folder  = '/media/veracrypt1/P03';
    calib1_file = 'KUKA_neckEMG_P03_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P03_AllOri2.txt';
elseif participant_num == 4
    rec_paths   = {'/media/veracrypt1/P04/2026-07-28_14-14-18/'};
    imu_folder  = '/media/veracrypt1/P04';
    calib1_file = 'KUKA_neckEMG_P04_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P04_AllOri2.txt';
elseif participant_num == 5
    rec_paths   = {'/media/veracrypt1/P05/2026-07-29_09-53-49/'};
    imu_folder  = '/media/veracrypt1/P05';
    calib1_file = 'KUKA_neckEMG_P05_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P05_AllOri2.txt';
elseif participant_num == 6
    rec_paths   = {'/media/veracrypt1/P06/2026-07-29_16-21-58/'};
    imu_folder  = '/media/veracrypt1/P06';
    calib1_file = 'KUKA_neckEMG_P06_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P06_AllOri2.txt';
elseif participant_num == 7
    rec_paths   = {'/media/veracrypt1/P07/2026-07-30_10-09-02/'};
    imu_folder  = '/media/veracrypt1/P07';
    calib1_file = 'KUKA_neckEMG_P07_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P07_AllOri2.txt';
elseif participant_num == 8
    rec_paths   = {'/media/veracrypt1/P08/2026-07-30_14-18-07/'};
    imu_folder  = '/media/veracrypt1/P08';
    calib1_file = 'KUKA_neckEMG_P08_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P08_AllOri2.txt';
elseif participant_num == 9
    rec_paths   = {'/media/veracrypt1/P09/2026-07-31_10-09-51/'};
    imu_folder  = '/media/veracrypt1/P09';
    calib1_file = 'KUKA_neckEMG_P09_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P09_AllOri2.txt';
elseif participant_num == 10
    rec_paths   = {'/media/veracrypt1/P10/2026-07-31_14-24-53/'};
    imu_folder  = '/media/veracrypt1/P10';
    calib1_file = 'KUKA_neckEMG_P10_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P10_AllOri2.txt';
end

% --- sync extractor's output folder for this participant (where it wrote
%     EMG_IMU_Match_Key.csv) -- must match its `outfile`/`csv_out_dir` ---
outfile     = sprintf('P%02d_sync_summary_order', participant_num);
csv_out_dir = fullfile(pwd, outfile);
key_path    = fullfile(csv_out_dir, 'EMG_IMU_Match_Key.csv');

% --- processing-specific settings _analysis_dataprep_aug22_SyncLowPassInterp_plot needs
%     (the sync extractor doesn't use these) ---
emg_sync_ch = 129;
emg_fs      = 2000;
imu_fs      = 4000;

% Every entry in muscle_list is run for every matched pair, so by default
% both electrode banks get processed (1 -> channels 1:64, 2 -> channels
% 65:128). Set to [1] or [2] to process only one bank.
muscle_list = [1, 2];

% --- direction is NO LONGER a fixed CONFIG value. It's read per-trial
%     from ClaudeKeySetting.m's IMU_Key_Table.csv (GroupLabel column),
%     joined onto the Match Key by IMU_FileName -- see LOAD THE KEY
%     TABLES section below. ---

% Folder where ClaudeKeySetting.m wrote IMU_Key_Table.csv for this
% participant. ClaudeKeySetting.m uses outfile = 'P%02d_sync_summary'
% (NOTE: no "_order" suffix, unlike csv_out_dir below for the Match Key --
% double check these are really two different folders on your setup,
% and adjust key_setting_out_dir here if you've since aligned the naming).
key_setting_outfile     = sprintf('P%02d_sync_summary', participant_num);
key_setting_out_dir     = fullfile(pwd, key_setting_outfile);
key_table_path          = fullfile(key_setting_out_dir, 'IMU_Key_Table.csv');

flange_calib_files = struct( ...
    'b_gyro', '/home/jslab/Documents/Matt/IMU Offsets/A1 (flange)/cal_20260810/b_gyr.txt', ...
    'b_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A1 (flange)/cal_20260810/b_acc.txt', ...
    's_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A1 (flange)/cal_20260810/S_acc.txt', ...
    't_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A1 (flange)/cal_20260810/T_acc.txt');

head_calib_files = struct( ...
    'b_gyro', '/home/jslab//Documents/Matt/IMU Offsets/A1 (flange)/cal_20260810/b_gyr.txt', ...
    'b_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A2 (head)/cal_20260810/b_acc.txt', ...
    's_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A2 (head)/cal_20260810/S_acc.txt', ...
    't_acc',  '/home/jslab/Documents/Matt/IMU Offsets/A2 (head)/cal_20260810/T_acc.txt');



%% ------------------------- LOAD THE MATCH KEY ---------------------------
if ~exist(key_path, 'file')
    error(['Match key not found: %s\n' ...
           'Run the sync pulse extractor for participant %d first.'], ...
           key_path, participant_num);
end

% Delimiter forced to ',' -- auto-detection can misfire here because
% IMU_FileName values are full of underscores (e.g.
% "KUKA_neckEMG_P05_Trial_Block1_01_02.txt"), and detectImportOptions can
% end up picking '_' as the delimiter instead of ',' since it looks more
% "consistent" across rows. Forcing it avoids that entirely.
key_table = readtable(key_path, 'TextType', 'string', 'Delimiter', ',');
n_pairs   = height(key_table);
fprintf('\nLoaded %d EMG-IMU pair(s) from:\n  %s\n', n_pairs, key_path);

required_key_cols = {'IMU_FileName', 'EMG_SessionIdx', 'EMG_LocalRec', ...
    'EMG_FirstHighToLowIdx', 'EMG_FinalLowToHighIdx', ...
    'IMU_FirstHighToLowIdx', 'IMU_FinalLowToHighIdx'};
missing_key_cols = setdiff(required_key_cols, key_table.Properties.VariableNames);
if ~isempty(missing_key_cols)
    error(['%s is missing expected column(s): %s\n' ...
           'Re-run the sync pulse extractor (it now writes the sync ' ...
           'indices directly into the Match Key).'], ...
           key_path, strjoin(missing_key_cols, ', '));
end

if n_pairs == 0
    fprintf('Nothing to process.\n');
    return;
end

%% ------------------------- LOAD THE DIRECTION KEY TABLE ------------------
% ClaudeKeySetting.m's IMU_Key_Table.csv: one row per IMU file, with a
% GroupLabel column ('ML' / 'AP' / 'DiaR' / 'DiaL') that maps 1:1 onto
% Claude_emg_imu_analysis_data_prep_aug18_v2's `direction` input. Join it onto the
% match key by filename so each row below knows its own direction.
if ~exist(key_table_path, 'file')
    error(['Direction key table not found: %s\n' ...
           'Run ClaudeKeySetting.m for participant %d first so every ' ...
           'IMU file has a GroupLabel (direction) assigned.'], ...
           key_table_path, participant_num);
end

% Delimiter forced to ',' for the same reason as key_path above --
% FileName values here are also full of underscores.
direction_key_table = readtable(key_table_path, 'TextType', 'string', 'Delimiter', ',');
required_cols = {'FileName', 'GroupLabel', 'ParameterSet', 'Condition'};
missing_cols  = setdiff(required_cols, direction_key_table.Properties.VariableNames);
if ~isempty(missing_cols)
    error('IMU_Key_Table.csv at %s is missing expected column(s): %s', ...
        key_table_path, strjoin(missing_cols, ', '));
end

% Map IMU_FileName -> GroupLabel / ParameterSet / Condition for fast per-row lookup.
fn_keys = cellstr(direction_key_table.FileName);
direction_map    = containers.Map(fn_keys, cellstr(direction_key_table.GroupLabel));
paramset_map     = containers.Map(fn_keys, cellstr(direction_key_table.ParameterSet));
condition_map    = containers.Map(fn_keys, cellstr(direction_key_table.Condition));

fprintf('Loaded direction key table (%d file(s)) from:\n  %s\n', ...
    height(direction_key_table), key_table_path);

n_muscles = numel(muscle_list);
n_runs    = n_pairs * n_muscles;
fprintf('Will run %d muscle setting(s) per pair -> %d total run(s).\n', ...
    n_muscles, n_runs);

%% ------------------------- BATCH LOOP -----------------------------------
log_ImuFile     = strings(n_runs, 1);
log_SessionIdx  = nan(n_runs, 1);
log_LocalRec    = nan(n_runs, 1);
log_Muscle      = nan(n_runs, 1);
log_Direction   = strings(n_runs, 1);
log_Status      = strings(n_runs, 1);
log_Message     = strings(n_runs, 1);
log_MatPath     = strings(n_runs, 1);
log_ElapsedSec  = nan(n_runs, 1);

run_i = 0;
for k = 1:n_pairs
    row         = key_table(k, :);   % this pair's full matched-trial row, incl. sync indices
    imu_file    = char(key_table.IMU_FileName(k));
    sess_idx    = key_table.EMG_SessionIdx(k);
    local_rec   = key_table.EMG_LocalRec(k);

    % --- Look up this trial's direction / ParameterSet / Condition from the key-setting table ---
    if isKey(direction_map, imu_file)
        this_direction  = string(direction_map(imu_file));
        this_paramset   = string(paramset_map(imu_file));
        this_condition  = string(condition_map(imu_file));
    else
        this_direction  = "";
        this_paramset   = "";
        this_condition  = "";
    end
    direction_missing = (this_direction == "" || ismissing(this_direction));

    for m = 1:n_muscles
        this_muscle = muscle_list(m);
        run_i = run_i + 1;

        fprintf('\n=== [%d/%d] Session %d, Rec %d  <->  %s  (Muscle %d, direction "%s", %s, %s) ===\n', ...
            run_i, n_runs, sess_idx, local_rec, imu_file, this_muscle, this_direction, this_paramset, this_condition);

        log_ImuFile(run_i)    = imu_file;
        log_SessionIdx(run_i) = sess_idx;
        log_LocalRec(run_i)   = local_rec;
        log_Muscle(run_i)     = this_muscle;
        log_Direction(run_i)  = this_direction;

        if sess_idx < 1 || sess_idx > numel(rec_paths)
            log_Status(run_i)  = "error";
            log_Message(run_i) = sprintf('EMG_SessionIdx %d out of range for rec_paths (numel=%d)', ...
                sess_idx, numel(rec_paths));
            fprintf('[ERROR] %s\n', log_Message(run_i));
            continue;
        end

        if direction_missing
            log_Status(run_i)  = "skipped_no_direction";
            log_Message(run_i) = sprintf('No GroupLabel found for "%s" in IMU_Key_Table.csv -- tag it with ClaudeKeySetting.m first.', imu_file);
            fprintf('[SKIPPED] %s\n', log_Message(run_i));
            continue;
        end

        P = struct();
        P.rec_path              = rec_paths{sess_idx};
        P.rec_num               = local_rec;
        P.imu_folder            = imu_folder;
        P.calib1_file           = calib1_file;
        P.calib2_file           = calib2_file;
        P.imu_file              = imu_file;
        P.emg_fs                = emg_fs;
        P.imu_fs                = imu_fs;
        P.emg_sync_ch           = emg_sync_ch;
        P.Muscle                = this_muscle;
        P.direction             = char(this_direction);
        P.ParameterSet          = char(this_paramset);
        P.condition             = char(this_condition);
        P.flange_calib_files    = flange_calib_files;
        P.head_calib_files      = head_calib_files;
        P.csv_out_dir           = csv_out_dir;

        % Added these with updated sync_extract_test Aug 22
        P.idxStart_emg          = row.EMG_FirstHighToLowIdx;
        P.idxEnd_emg            = row.EMG_FinalLowToHighIdx;
        P.idxStart_imu          = row.IMU_FirstHighToLowIdx;
        P.idxEnd_imu            = row.IMU_FinalLowToHighIdx;

        P.make_diagnostic_plots = true;
        % P.expected_emg/P.expected_imu below are NOT available here --
        % EMG_Pre/EMG_Trial/EMG_Post/IMU_Pre/IMU_Trial/IMU_Post only exist
        % in EMG_IMU_Matched_Trials.csv, not in the (deliberately trimmer)
        % Match Key. Load that file separately (joined on IMU_FileName) if
        % you want these re-enabled.
        % P.expected_emg = struct('PreMarkerDur', row.EMG_Pre, 'TrialDur', row.EMG_Trial, 'PostMarkerDur', row.EMG_Post);
        % P.expected_imu = struct('PreMarkerDur', row.IMU_Pre, 'TrialDur', row.IMU_Trial, 'PostMarkerDur', row.IMU_Post);

        % NOTE: calls _analysis_dataprep_aug22_SyncLowPassInterp_plot -- must match
        % the actual function name in that .m file (MATLAB requires the
        % function name and filename to match). If you rename that file,
        % update this call to match.
        t_start = tic;
        try
            data = KUKA_analysis_preprocess_consumer(P);
            log_ElapsedSec(run_i) = toc(t_start);

            if isfield(data, 'status') && strcmp(data.status, 'skipped_sync_fail')
                log_Status(run_i)  = "skipped_sync_fail";
                log_Message(run_i) = "Sync validation failed (see Skipped_Sync_Trials.csv)";
                fprintf('[SKIPPED] %s\n', log_Message(run_i));
            else
                log_Status(run_i) = "ok";
                if isfield(data, 'meta') && isfield(data.meta, 'mat_path')
                    log_MatPath(run_i) = data.meta.mat_path;
                end
                fprintf('[OK] %.1f s\n', log_ElapsedSec(run_i));
            end
        catch ME
            log_ElapsedSec(run_i) = toc(t_start);
            log_Status(run_i)  = "error";
            log_Message(run_i) = string(ME.message);
            fprintf('[ERROR] %s\n', ME.message);
        end
    end
end

%% ------------------------- SUMMARY & LOG --------------------------------
batch_log = table(log_ImuFile, log_SessionIdx, log_LocalRec, log_Muscle, ...
    log_Direction, log_Status, log_Message, log_MatPath, log_ElapsedSec, ...
    'VariableNames', {'IMU_FileName', 'EMG_SessionIdx', 'EMG_LocalRec', ...
    'Muscle', 'Direction', 'Status', 'Message', 'MatPath', 'ElapsedSec'});

disp(batch_log);

log_out_path = fullfile(csv_out_dir, 'Batch_Data_Prep_Log.csv');
writetable(batch_log, log_out_path);

n_ok               = sum(log_Status == "ok");
n_skipped_sync      = sum(log_Status == "skipped_sync_fail");
n_skipped_direction = sum(log_Status == "skipped_no_direction");
n_error             = sum(log_Status == "error");
fprintf(['\n=== BATCH COMPLETE: %d ok, %d skipped (sync fail), ' ...
    '%d skipped (no direction), %d error, out of %d run(s) [%d pair(s) x %d muscle setting(s)] ===\n'], ...
    n_ok, n_skipped_sync, n_skipped_direction, n_error, n_runs, n_pairs, n_muscles);
fprintf('Log written to: %s\n', log_out_path);