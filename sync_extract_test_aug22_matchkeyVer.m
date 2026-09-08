%% ========================================================================
%  EMG & IMU SYNC PULSE EXTRACTOR AND TRIAL MATCHER  (v6.1)
% ========================================================================
%  WHAT CHANGED VS. v6
%  --------------------
%  1. IMU Index Extraction: Tracked and stored IMU_FirstHighToLowIdx 
%     and IMU_FinalLowToHighIdx in the output tables.
%  2. Match Output Structs: Updated candidate assignment and matching 
%     functions (ordered and unordered) to propagate IMU sample indices.
% ========================================================================

clear; clc; close all;

%% ------------------------- CONFIG -------------------------------------
participant_num = 10; % Set to the participant you are running

if participant_num == 1
    rec_paths = {'/media/veracrypt1/P01/2026-07-24_14-02-02/'};
    imu_folder = '/media/veracrypt1/P01';
    calib1_file = 'KUKA_neckEMG_P01_Trial_calib1_redo2.txt';
    calib2_file = 'KUKA_neckEMG_P01_Trial_calib2.txt';
elseif participant_num == 2
    % Multiple sessions for P02
    rec_paths = {
        '/media/veracrypt1/P02/2026-07-27_10-15-17/', ...
        '/media/veracrypt1/P02/2026-07-27_12-03-56/'
    };
    imu_folder = '/media/veracrypt1/P02'; 
    calib1_file = 'KUKA_neckEMG_P02_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P02_AllOri2.txt';
elseif participant_num == 3
    rec_paths = {'/media/veracrypt1/P03/2026-07-27_14-37-14/'};
    imu_folder = '/media/veracrypt1/P03'; 
    calib1_file = 'KUKA_neckEMG_P03_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P03_AllOri2.txt';
elseif participant_num == 4
    rec_paths = {'/media/veracrypt1/P04/2026-07-28_14-14-18/'};
    imu_folder = '/media/veracrypt1/P04'; 
    calib1_file = 'KUKA_neckEMG_P04_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P04_AllOri2.txt';
elseif participant_num == 5
    rec_paths = {'/media/veracrypt1/P05/2026-07-29_09-53-49/'};
    imu_folder = '/media/veracrypt1/P05';
    calib1_file = 'KUKA_neckEMG_P05_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P05_AllOri2.txt';
elseif participant_num == 6
    rec_paths = {'/media/veracrypt1/P06/2026-07-29_16-21-58/'};
    imu_folder = '/media/veracrypt1/P06';
    calib1_file = 'KUKA_neckEMG_P06_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P06_AllOri2.txt';
elseif participant_num == 7
    rec_paths = {'/media/veracrypt1/P07/2026-07-30_10-09-02/'};
    imu_folder = '/media/veracrypt1/P07';
    calib1_file = 'KUKA_neckEMG_P07_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P07_AllOri2.txt';
elseif participant_num == 8
    rec_paths = {'/media/veracrypt1/P08/2026-07-30_14-18-07/'};
    imu_folder = '/media/veracrypt1/P08';
    calib1_file = 'KUKA_neckEMG_P08_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P08_AllOri2.txt';
elseif participant_num == 9
    rec_paths = {'/media/veracrypt1/P09/2026-07-31_10-09-51/'};
    imu_folder = '/media/veracrypt1/P09';
    calib1_file = 'KUKA_neckEMG_P09_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P09_AllOri2.txt';
elseif participant_num == 10
    rec_paths = {'/media/veracrypt1/P10/2026-07-31_14-24-53/'};
    imu_folder = '/media/veracrypt1/P10';
    calib1_file = 'KUKA_neckEMG_P10_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P10_AllOri2.txt';
end

emg_sync_ch = 129;      % Sync channel index in OpenEphys stream
emg_fs      = 2000;     % EMG sampling frequency (Hz)
imu_fs      = 4000;     % IMU sampling frequency (Hz)

cfg = struct();
cfg.hyst_low_frac    = 0.40;   
cfg.hyst_high_frac   = 0.60;   
cfg.min_glitch_time  = 0.005;  
cfg.min_trial_dur    = 15.0;   
cfg.max_marker_dur   = 15.0;   
cfg.match_tol        = 0.05; 
cfg.min_high_dur     = 0.5;
cfg.match_mode       = 'unordered'; 

outfile = sprintf('P%02d_sync_summary_order', participant_num);
csv_out_dir = fullfile(pwd, outfile);
qc_plot_dir = fullfile(csv_out_dir, 'QC_plots');
if ~exist(csv_out_dir, 'dir'); mkdir(csv_out_dir); end
if ~exist(qc_plot_dir, 'dir'); mkdir(qc_plot_dir); end

make_qc_plots = false; 

%% 1. Process OpenEphys EMG recordings ==================================
emg_results  = struct([]);
arr_idx = 0; % Internal continuous index just for array sizing

for s = 1:numel(rec_paths)
    rec_path = rec_paths{s};
    fprintf('\n=== Processing EMG Session %d / %d : %s ===\n', s, numel(rec_paths), rec_path);
    
    session = Session(char(rec_path));
    node    = session.recordNodes{1};
    num_emg_recs = length(node.recordings);

    for r = 1:num_emg_recs
        rec = node.recordings{r};
        streams = rec.continuous.keys();
        if isempty(streams)
            continue;
        end
        
        arr_idx = arr_idx + 1; % Increment only for valid recordings
        
        stream_name = streams{1};
        stream = rec.continuous(stream_name);

        exp_num = 1;
        try
            if isprop(rec, 'experimentNumber') || isfield(rec, 'experimentNumber')
                exp_num = rec.experimentNumber;
            end
        catch
        end

        node_name = 'Node_1';
        try
            if isprop(node, 'name') && ~isempty(node.name)
                node_name = node.name;
            end
        catch
        end

        bit_volts = [];
        try
            if iscell(stream.channels)
                bit_volts = stream.channels{emg_sync_ch}.bit_volts;
            elseif isstruct(stream.channels)
                bit_volts = stream.channels(emg_sync_ch).bit_volts;
            end
        catch
        end
        if isempty(bit_volts) && isfield(stream, 'header') && isfield(stream.header, 'bitVolts')
            bit_volts = stream.header.bitVolts;
        end
        if isempty(bit_volts)
            bit_volts = 0.195; % uV/bit default Intan factor
        end

        raw_samples  = double(stream.samples(emg_sync_ch, :));
        sync_sig_mV  = (raw_samples * bit_volts) / 1000;
        t_vec        = (0:length(sync_sig_mV)-1) / emg_fs;

        % Name reflects Session and Local Recording
        rec_name = sprintf('Sess%d_Rec%d', s, r); 
        segments = extract_all_segments(t_vec, sync_sig_mV, emg_fs, cfg);
        trials   = identify_trials(segments, cfg);

        if make_qc_plots
            plot_sync_qc(t_vec, sync_sig_mV, segments, trials, sprintf('EMG %s', rec_name), ...
                fullfile(qc_plot_dir, sprintf('EMG_%s.png', rec_name)));
        end

        % Store the two-part key
        emg_results(arr_idx).SessionIdx     = s;
        emg_results(arr_idx).LocalRecNumber = r;
        emg_results(arr_idx).RecName        = rec_name;
        emg_results(arr_idx).Node           = node_name;
        emg_results(arr_idx).Stream         = char(stream_name);
        emg_results(arr_idx).ExpNumber      = exp_num;
        emg_results(arr_idx).FullPath       = fullfile(rec_path, sprintf('recording_%d', r));
        emg_results(arr_idx).TotalSamples   = length(sync_sig_mV);
        emg_results(arr_idx).DurationSec    = t_vec(end);
        emg_results(arr_idx).Segments       = segments;
        emg_results(arr_idx).Trials         = trials;
        emg_results(arr_idx).NumTrials      = numel(trials);

        if numel(trials) > 1
            fprintf('NOTE: %s has %d confirmed low-high-low trial triplets - expected 1. Inspect its QC plot.\n', ...
                rec_name, numel(trials));
        elseif numel(trials) == 0
            fprintf('NOTE: %s has NO confirmed trial triplet. Check its QC plot.\n', rec_name);
        end
    end
end

%% 2. Process text-based IMU trial files =================================
imu_files = dir(fullfile(imu_folder, 'KUKA_*.txt'));
num_imu_files = length(imu_files);

imu_results = struct([]);
skipped_imu_results = struct([]);
valid_count = 0; skipped_count = 0;

for f = 1:num_imu_files
    fname = imu_files(f).name;
    fpath = fullfile(imu_folder, fname);

    try
        tbl = load_imu_table_v2(fpath);
    catch ME
        fprintf('--- ERROR IN FILE: %s ---\n', fname);
        disp(ME.getReport());
        tbl = table();
    end
    [nRows, nCols] = size(tbl);

    if nRows < 10 || nCols < 12
        skipped_count = skipped_count + 1;
        skipped_imu_results(skipped_count).FileName = fname;
        skipped_imu_results(skipped_count).Rows     = nRows;
        skipped_imu_results(skipped_count).Cols     = nCols;
        skipped_imu_results(skipped_count).Reason   = sprintf('Dimensional threshold missed (%dx%d)', nRows, nCols);
        continue;
    end

    valid_count = valid_count + 1;

    if ismember('Sync', tbl.Properties.VariableNames)
        sync_sig = tbl.Sync;
    else
        sync_arr = table2array(tbl);
        sync_sig = sync_arr(:, end);
    end

    if ismember('Time', tbl.Properties.VariableNames)
        if iscell(tbl.Time)
            t_vec = str2double(tbl.Time);
        else
            t_vec = tbl.Time;
        end
    else
        t_vec = (0:length(sync_sig)-1) / imu_fs;
    end
    t_vec = t_vec - t_vec(1);

    segments = extract_all_segments(t_vec, sync_sig, imu_fs, cfg);
    trials   = identify_trials(segments, cfg);

    if make_qc_plots
        [~, base, ~] = fileparts(fname);
        plot_sync_qc(t_vec, sync_sig, segments, trials, sprintf('IMU %s', base), ...
            fullfile(qc_plot_dir, sprintf('IMU_%s.png', base)));
    end

    imu_results(valid_count).FileName      = fname;
    imu_results(valid_count).TotalSamples  = length(sync_sig);
    imu_results(valid_count).DurationSec   = t_vec(end);
    imu_results(valid_count).Segments      = segments;
    imu_results(valid_count).Trials        = trials;
    imu_results(valid_count).NumTrials     = numel(trials);

    if numel(trials) > 1
        fprintf('NOTE: %s has %d confirmed low-high-low trial triplets - expected 1 per file. Inspect its QC plot.\n', ...
            fname, numel(trials));
    elseif numel(trials) == 0
        fprintf('NOTE: %s has NO confirmed trial triplet. Check its QC plot.\n', fname);
    end
end

%% 3. Build & export per-file/recording summary tables ===================
fprintf('\n=== EMG RECORDINGS SYNC SUMMARY ===\n');
EMG_SessionIdx = [emg_results.SessionIdx]';
EMG_LocalRec   = [emg_results.LocalRecNumber]';
EMG_NumTrials  = [emg_results.NumTrials]';
emg_table = table(EMG_SessionIdx, EMG_LocalRec, EMG_NumTrials);
disp(emg_table);
writetable(emg_table, fullfile(csv_out_dir, 'EMG_Sync_Summary.csv'));

fprintf('\n=== IMU FILES SYNC SUMMARY ===\n');
if valid_count > 0
    IMU_FileName  = {imu_results.FileName}';
    IMU_NumTrials = [imu_results.NumTrials]';
    imu_table = table(IMU_FileName, IMU_NumTrials);
    disp(imu_table);
    writetable(imu_table, fullfile(csv_out_dir, 'IMU_Sync_Summary.csv'));
else
    disp('No valid IMU files met the row/column size criteria.');
end

if skipped_count > 0
    Skipped_FileName = {skipped_imu_results.FileName}';
    Skipped_Rows = [skipped_imu_results.Rows]';
    Skipped_Cols = [skipped_imu_results.Cols]';
    Skipped_Reason = {skipped_imu_results.Reason}';
    skipped_table = table(Skipped_FileName, Skipped_Rows, Skipped_Cols, Skipped_Reason);
    writetable(skipped_table, fullfile(csv_out_dir, 'Skipped_IMU_Files.csv'));
end

imu_block_all = nan(numel(imu_results), 1);
imu_trial_all = nan(numel(imu_results), 1);
for i = 1:numel(imu_results)
    [imu_block_all(i), imu_trial_all(i)] = parse_block_trial(imu_results(i).FileName);
end
[~, file_sort_ix] = sortrows([imu_block_all, imu_trial_all], [1 2]);
for rank = 1:numel(file_sort_ix)
    imu_results(file_sort_ix(rank)).FileIndex = rank;
end

%% 4. Match EMG trials <-> IMU trials =====================================
fprintf('\n=== EMG-IMU TRIAL MATCHING (mode=%s, full triplet, tol=%.3fs/component) ===\n', cfg.match_mode, cfg.match_tol);

emg_cand = struct('EmgIdx', {}, 'SessionIdx', {}, 'LocalRecNumber', {}, 'TrialIdxInRec', {}, ...
    'NumTrialsInRec', {}, 'PreMarkerDur', {}, 'TrialDur', {}, 'PostMarkerDur', {}, ...
    'PreStartIdx', {}, 'PostEndIdx', {}, 'ChronoIdx', {});
for e = 1:numel(emg_results)
    for k = 1:emg_results(e).NumTrials
        tr = emg_results(e).Trials(k);
        emg_cand(end+1) = struct('EmgIdx', e, ...
            'SessionIdx', emg_results(e).SessionIdx, ...
            'LocalRecNumber', emg_results(e).LocalRecNumber, ...
            'TrialIdxInRec', k, ...
            'NumTrialsInRec', emg_results(e).NumTrials, ... % Tracked for conditional matching
            'PreMarkerDur', tr.PreMarkerDur, 'TrialDur', tr.TrialDur, 'PostMarkerDur', tr.PostMarkerDur, ...
            'PreStartIdx', tr.PreStartIdx, 'PostEndIdx', tr.PostEndIdx, ...
            'ChronoIdx', numel(emg_cand)+1); %#ok<SAGROW>
    end
end

imu_cand = struct('ImuIdx', {}, 'TrialIdxInFile', {}, 'PreMarkerDur', {}, 'TrialDur', {}, 'PostMarkerDur', {}, ...
    'PreStartIdx', {}, 'PostEndIdx', {}, 'BlockNum', {}, 'TrialNum', {}, 'ChronoIdx', {}, 'FileIndex', {});
n_unparsed = 0;
for i = 1:numel(imu_results)
    [blk, trl] = parse_block_trial(imu_results(i).FileName);
    if isnan(blk)
        n_unparsed = n_unparsed + 1;
    end
    for k = 1:imu_results(i).NumTrials
        tr = imu_results(i).Trials(k);
        imu_cand(end+1) = struct('ImuIdx', i, 'TrialIdxInFile', k, ...
            'PreMarkerDur', tr.PreMarkerDur, 'TrialDur', tr.TrialDur, 'PostMarkerDur', tr.PostMarkerDur, ...
            'PreStartIdx', tr.PreStartIdx, 'PostEndIdx', tr.PostEndIdx, ...
            'BlockNum', blk, 'TrialNum', trl, 'ChronoIdx', NaN, 'FileIndex', imu_results(i).FileIndex); %#ok<SAGROW>
    end
end
if n_unparsed > 0
    fprintf('WARNING: could not parse block/trial number from %d IMU filename(s) - they were sorted last.\n', n_unparsed);
end

if ~isempty(imu_cand)
    [~, sort_ix] = sortrows([[imu_cand.BlockNum]', [imu_cand.TrialNum]'], [1 2]);
    imu_cand = imu_cand(sort_ix);
    for c = 1:numel(imu_cand)
        imu_cand(c).ChronoIdx = c;
    end
end

switch cfg.match_mode
    case 'ordered'
        [match_results, unmatched_emg, unmatched_imu] = match_trials_ordered(emg_cand, imu_cand, cfg.match_tol);
    case 'unordered'
        [match_results, unmatched_emg, unmatched_imu] = match_trials_unordered(emg_cand, imu_cand, cfg.match_tol);
    otherwise
        error('cfg.match_mode must be ''ordered'' or ''unordered''');
end

if ~isempty(match_results)
    [~, sort_by_emg] = sort([match_results.EmgChronoIdx]);
    match_results = match_results(sort_by_emg);
    imu_chrono_seq = [match_results.ImuChronoIdx];
    in_order_mask = compute_lis_mask(imu_chrono_seq);
    for m = 1:numel(match_results)
        match_results(m).InOrder = in_order_mask(m);
    end
end

if ~isempty(match_results)
    n = numel(match_results);
    EMG_SessionIdx = zeros(n,1); EMG_LocalRec = zeros(n,1); 
    EMG_FullPath = cell(n,1); IMU_FileName = cell(n,1);
    EMG_ChronoIdx = zeros(n,1); IMU_ChronoIdx = zeros(n,1); IMU_FileIndex = zeros(n,1);
    EMG_FirstHighToLowIdx = zeros(n,1); EMG_FinalLowToHighIdx = zeros(n,1);
    IMU_FirstHighToLowIdx = zeros(n,1); IMU_FinalLowToHighIdx = zeros(n,1);
    EMG_Pre = zeros(n,1); EMG_Trial = zeros(n,1); EMG_Post = zeros(n,1);
    IMU_Pre = zeros(n,1); IMU_Trial = zeros(n,1); IMU_Post = zeros(n,1);
    PreDiff_ms = zeros(n,1); TrialDiff_ms = zeros(n,1); PostDiff_ms = zeros(n,1);
    InOrder = false(n,1);
    
    for m = 1:n
        e = match_results(m).EmgIdx; i = match_results(m).ImuIdx;
        
        EMG_SessionIdx(m) = emg_results(e).SessionIdx;
        EMG_LocalRec(m)   = emg_results(e).LocalRecNumber;
        
        EMG_FullPath{m}  = emg_results(e).FullPath;
        IMU_FileName{m}  = imu_results(i).FileName;
        EMG_ChronoIdx(m) = match_results(m).EmgChronoIdx;
        IMU_ChronoIdx(m) = match_results(m).ImuChronoIdx;
        IMU_FileIndex(m) = match_results(m).ImuFileIndex;
        
        % Extract exact EMG sample indices
        EMG_FirstHighToLowIdx(m) = match_results(m).EmgPreStartIdx;
        EMG_FinalLowToHighIdx(m) = match_results(m).EmgPostEndIdx;
        
        % Extract exact IMU sample indices
        IMU_FirstHighToLowIdx(m) = match_results(m).ImuPreStartIdx;
        IMU_FinalLowToHighIdx(m) = match_results(m).ImuPostEndIdx;
        
        EMG_Pre(m) = match_results(m).EmgPre; EMG_Trial(m) = match_results(m).EmgTrial; EMG_Post(m) = match_results(m).EmgPost;
        IMU_Pre(m) = match_results(m).ImuPre; IMU_Trial(m) = match_results(m).ImuTrial; IMU_Post(m) = match_results(m).ImuPost;
        PreDiff_ms(m)   = (IMU_Pre(m)   - EMG_Pre(m))   * 1000;
        TrialDiff_ms(m) = (IMU_Trial(m) - EMG_Trial(m)) * 1000;
        PostDiff_ms(m)  = (IMU_Post(m)  - EMG_Post(m))  * 1000;
        InOrder(m) = match_results(m).InOrder;
    end
    
    match_table = table(EMG_SessionIdx, EMG_LocalRec, EMG_ChronoIdx, IMU_ChronoIdx, IMU_FileIndex, ...
        EMG_FirstHighToLowIdx, EMG_FinalLowToHighIdx, ...
        IMU_FirstHighToLowIdx, IMU_FinalLowToHighIdx, ...
        EMG_FullPath, IMU_FileName, ...
        EMG_Pre, EMG_Trial, EMG_Post, IMU_Pre, IMU_Trial, IMU_Post, ...
        PreDiff_ms, TrialDiff_ms, PostDiff_ms, InOrder);
    match_table = sortrows(match_table, {'EMG_SessionIdx', 'EMG_LocalRec'});
    disp(match_table);
    writetable(match_table, fullfile(csv_out_dir, 'EMG_IMU_Matched_Trials.csv'));

    n_out_of_order = sum(~InOrder);
    fprintf('\n%d of %d matched trials are NOT consistent with a single chronological ordering (InOrder=false).\n', ...
        n_out_of_order, n);

    % Clean matching key for downstream use -- also carries the sync
    % sample indices (EMG_FirstHighToLowIdx/EMG_FinalLowToHighIdx/
    % IMU_FirstHighToLowIdx/IMU_FinalLowToHighIdx) so downstream analysis
    % code can read this one file instead of separately loading
    % EMG_IMU_Matched_Trials.csv.
    match_key = table(IMU_FileName, IMU_FileIndex, EMG_SessionIdx, EMG_LocalRec, EMG_FullPath, ...
        EMG_FirstHighToLowIdx, EMG_FinalLowToHighIdx, ...
        IMU_FirstHighToLowIdx, IMU_FinalLowToHighIdx);
    match_key = sortrows(match_key, 'IMU_FileIndex');

    key_table_path = fullfile(csv_out_dir, 'IMU_Key_Table.csv');
    if exist(key_table_path, 'file')
        param_key = readtable(key_table_path, 'TextType', 'string');
        match_key.IMU_FileName = string(match_key.IMU_FileName);
        param_key.FileName = string(param_key.FileName);
        if all(ismember({'ParameterSet','Condition'}, param_key.Properties.VariableNames))
            match_key = outerjoin(match_key, param_key(:, {'FileName','ParameterSet','Condition'}), ...
                'LeftKeys', 'IMU_FileName', 'RightKeys', 'FileName', 'MergeKeys', true, 'Type', 'left');
            fprintf('\nMerged ParameterSet/Condition from %s into the matching key.\n', key_table_path);
        end
    else
        fprintf('\nNOTE: %s not found - matching key has no ParameterSet/Condition columns.\n', key_table_path);
    end

    match_key_path = fullfile(csv_out_dir, 'EMG_IMU_Match_Key.csv');
    writetable(match_key, match_key_path);
    fprintf('\nClean matching key saved to: %s\n', match_key_path);
else
    fprintf('No matches within tolerance.\n');
end

if ~isempty(unmatched_imu)
    fprintf('\n--- UNMATCHED IMU TRIALS (closest EMG candidate by trial duration, for diagnosis) ---\n');
    n = numel(unmatched_imu);
    IMU_FileName = cell(n,1); IMU_Pre = zeros(n,1); IMU_Trial = zeros(n,1); IMU_Post = zeros(n,1);
    Closest_EMG_Session = zeros(n,1); Closest_EMG_LocalRec = zeros(n,1); 
    Closest_EMG_Pre = zeros(n,1); Closest_EMG_Trial = zeros(n,1); Closest_EMG_Post = zeros(n,1);
    
    for u = 1:n
        i = unmatched_imu(u).ImuIdx;
        IMU_FileName{u} = imu_results(i).FileName;
        IMU_Pre(u) = unmatched_imu(u).PreMarkerDur;
        IMU_Trial(u) = unmatched_imu(u).TrialDur;
        IMU_Post(u) = unmatched_imu(u).PostMarkerDur;
        if ~isempty(emg_cand)
            [~, jmin] = min(abs([emg_cand.TrialDur] - unmatched_imu(u).TrialDur));
            Closest_EMG_Session(u) = emg_results(emg_cand(jmin).EmgIdx).SessionIdx;
            Closest_EMG_LocalRec(u) = emg_results(emg_cand(jmin).EmgIdx).LocalRecNumber;
            Closest_EMG_Pre(u)   = emg_cand(jmin).PreMarkerDur;
            Closest_EMG_Trial(u) = emg_cand(jmin).TrialDur;
            Closest_EMG_Post(u)  = emg_cand(jmin).PostMarkerDur;
        else
            Closest_EMG_Session(u) = NaN; Closest_EMG_LocalRec(u) = NaN; 
            Closest_EMG_Pre(u) = NaN; Closest_EMG_Trial(u) = NaN; Closest_EMG_Post(u) = NaN;
        end
    end
    unmatched_table = table(IMU_FileName, IMU_Pre, IMU_Trial, IMU_Post, ...
        Closest_EMG_Session, Closest_EMG_LocalRec, Closest_EMG_Pre, Closest_EMG_Trial, Closest_EMG_Post);
    disp(unmatched_table);
    writetable(unmatched_table, fullfile(csv_out_dir, 'Unmatched_IMU_Trials.csv'));
else
    fprintf('\nAll detected IMU trials matched an EMG trial within tolerance.\n');
end

if ~isempty(unmatched_emg)
    fprintf('\n--- UNMATCHED EMG TRIALS (%d) ---\n', numel(unmatched_emg));
    for u = 1:numel(unmatched_emg)
        e = unmatched_emg(u).EmgIdx;
        fprintf('  Session %d, Recording %d (%s): trial %.3fs (pre %.3fs / post %.3fs) had no IMU match within tolerance\n', ...
            emg_results(e).SessionIdx, emg_results(e).LocalRecNumber, emg_results(e).FullPath, ...
            unmatched_emg(u).TrialDur, unmatched_emg(u).PreMarkerDur, unmatched_emg(u).PostMarkerDur);
    end
end

fprintf('\nQC plots (one per recording/file) saved to: %s\n', qc_plot_dir);


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function segments = extract_all_segments(t_vec, sig, fs, cfg)
    sig = double(sig(:));
    t_vec = double(t_vec(:));
    n = numel(sig);

    segments = struct('Type', {}, 'StartTime', {}, 'EndTime', {}, 'Duration', {}, 'StartIdx', {}, 'EndIdx', {});

    % 1. Determine hysteresis thresholds
    sorted_sig = sort(sig);
    k = max(1, round(0.10 * n));
    val_low  = mean(sorted_sig(1:k));
    val_high = mean(sorted_sig(end-k+1:end));
    span = val_high - val_low;
    if span <= 0 || ~isfinite(span)
        return;
    end
    th_low  = val_low + cfg.hyst_low_frac  * span;
    th_high = val_low + cfg.hyst_high_frac * span;

    % 2. Calculate initial state transitions
    above = sig > th_high;
    below = sig < th_low;
    raw = zeros(n, 1);
    raw(above) =  1;
    raw(below) = -1;
    cross_idx = find(raw ~= 0);
    init_state = sig(1) > th_high;
    if isempty(cross_idx)
        state = repmat(init_state, n, 1);
    else
        fillidx = zeros(n, 1);
        fillidx(cross_idx) = cross_idx;
        fillidx = cummax(fillidx);
        has_crossed = fillidx > 0;
        if cross_idx(1) == 1
            has_crossed(1) = true;
        end
        state = false(n, 1);
        state(has_crossed)  = raw(fillidx(has_crossed)) > 0;
        state(~has_crossed) = init_state;
    end

    % 3. Populate raw segments struct (Preserving i0 and i1 boundary sample indices)
    change_idx = find(diff(state) ~= 0);
    bounds = [1; change_idx + 1; n + 1]; 
    n_seg = numel(bounds) - 1;
    for s = 1:n_seg
        i0 = bounds(s); i1 = bounds(s+1) - 1;
        if state(i0)
            typ = 'high';
        else
            typ = 'low';
        end
        segments(end+1) = struct('Type', typ, 'StartTime', t_vec(i0), 'EndTime', t_vec(i1), ...
            'Duration', t_vec(i1) - t_vec(i0), 'StartIdx', i0, 'EndIdx', i1); %#ok<AGROW>
    end

    % 4. Filter glitch durations and false high states
    if isfield(cfg, 'min_high_dur')
        min_high_dur = cfg.min_high_dur;
    end

    changed = true;
    while changed && numel(segments) > 1
        changed = false;
        
        % Clean up sub-glitch transients (updating sample index bounds on merges)
        durations = [segments.Duration];
        tiny_ix = find(durations < cfg.min_glitch_time, 1, 'first');
        
        if ~isempty(tiny_ix)
            if tiny_ix == 1
                segments(2).StartTime = segments(1).StartTime;
                segments(2).StartIdx  = segments(1).StartIdx;
                segments(2).Duration  = segments(2).EndTime - segments(2).StartTime;
                segments(1) = [];
            elseif tiny_ix == numel(segments)
                segments(end-1).EndTime  = segments(end).EndTime;
                segments(end-1).EndIdx   = segments(end).EndIdx;
                segments(end-1).Duration = segments(end-1).EndTime - segments(end-1).StartTime;
                segments(end) = [];
            else
                segments(tiny_ix-1).EndTime  = segments(tiny_ix+1).EndTime;
                segments(tiny_ix-1).EndIdx   = segments(tiny_ix+1).EndIdx;
                segments(tiny_ix-1).Duration = segments(tiny_ix-1).EndTime - segments(tiny_ix-1).StartTime;
                segments([tiny_ix, tiny_ix+1]) = [];
            end
            changed = true;
            continue; 
        end

        % Filter out false 'high' transitions < min_high_dur
        for i = 1:numel(segments)
            if strcmp(segments(i).Type, 'high') && segments(i).Duration < min_high_dur
                if i > 1 && i < numel(segments)
                    segments(i-1).EndTime  = segments(i+1).EndTime;
                    segments(i-1).EndIdx   = segments(i+1).EndIdx;
                    segments(i-1).Duration = segments(i-1).EndTime - segments(i-1).StartTime;
                    segments([i, i+1]) = [];
                elseif i == 1 && numel(segments) > 1
                    segments(1).Type = 'low';
                elseif i == numel(segments) && numel(segments) > 1
                    segments(i-1).EndTime  = segments(i).EndTime;
                    segments(i-1).EndIdx   = segments(i).EndIdx;
                    segments(i-1).Duration = segments(i-1).EndTime - segments(i-1).StartTime;
                    segments(i) = [];
                end
                changed = true;
                break; 
            end
        end
    end
end

function trials = identify_trials(segments, cfg)
    trials = struct('TrialIdx', {}, 'PreMarkerDur', {}, 'TrialDur', {}, 'PostMarkerDur', {}, ...
        'TrialStart', {}, 'TrialEnd', {}, 'PreMarkerStart', {}, 'PostMarkerEnd', {}, ...
        'PreStartIdx', {}, 'PostEndIdx', {});
    n = numel(segments);
    for k = 1:n
        if ~strcmp(segments(k).Type, 'high') || segments(k).Duration < cfg.min_trial_dur
            continue;
        end
        has_pre  = (k > 1) && strcmp(segments(k-1).Type, 'low') && segments(k-1).Duration <= cfg.max_marker_dur;
        has_post = (k < n) && strcmp(segments(k+1).Type, 'low') && segments(k+1).Duration <= cfg.max_marker_dur;
        if has_pre && has_post
            trials(end+1) = struct('TrialIdx', k, ...
                'PreMarkerDur', segments(k-1).Duration, 'TrialDur', segments(k).Duration, 'PostMarkerDur', segments(k+1).Duration, ...
                'TrialStart', segments(k).StartTime, 'TrialEnd', segments(k).EndTime, ...
                'PreMarkerStart', segments(k-1).StartTime, 'PostMarkerEnd', segments(k+1).EndTime, ...
                'PreStartIdx', segments(k-1).StartIdx, ... % Direct index of first high-to-low
                'PostEndIdx', segments(k+1).EndIdx);       % Direct index of final low-to-high %#ok<AGROW>
        end
    end
end


function [blockNum, trialNum] = parse_block_trial(fname)
    fname = char(fname);
    tok = regexpi(fname, 'block(\d+)_(\d+)', 'tokens', 'once');
    if isempty(tok)
        blockNum = NaN; trialNum = NaN;
    else
        blockNum = str2double(tok{1});
        trialNum = str2double(tok{2});
    end
end


function [matches, unmatched_emg, unmatched_imu] = match_trials_ordered(emg_cand, imu_cand, tol)
    n = numel(emg_cand);
    m = numel(imu_cand);
    matches = struct('EmgIdx', {}, 'ImuIdx', {}, 'EmgPre', {}, 'EmgTrial', {}, 'EmgPost', {}, ...
        'ImuPre', {}, 'ImuTrial', {}, 'ImuPost', {}, 'EmgPreStartIdx', {}, 'EmgPostEndIdx', {}, ...
        'ImuPreStartIdx', {}, 'ImuPostEndIdx', {}, ...
        'EmgChronoIdx', {}, 'ImuChronoIdx', {}, 'ImuFileIndex', {}, 'InOrder', {});

    if n == 0 || m == 0
        unmatched_emg = emg_cand;
        unmatched_imu = imu_cand;
        return;
    end

    skip_cost = 3 * tol;  
    big = 1e9;

    C = zeros(n+1, m+1);
    ptr = zeros(n+1, m+1); 
    for i = 1:n
        C(i+1,1) = C(i,1) + skip_cost;
        ptr(i+1,1) = 2;
    end
    for j = 1:m
        C(1,j+1) = C(1,j) + skip_cost;
        ptr(1,j+1) = 3;
    end

    for i = 1:n
        for j = 1:m
            d_pre   = abs(emg_cand(i).PreMarkerDur  - imu_cand(j).PreMarkerDur);
            d_trial = abs(emg_cand(i).TrialDur       - imu_cand(j).TrialDur);
            d_post  = abs(emg_cand(i).PostMarkerDur - imu_cand(j).PostMarkerDur);
            if d_pre <= tol && d_trial <= tol && d_post <= tol
                match_cost = d_pre + d_trial + d_post;
            else
                match_cost = big;
            end
            diag_c = C(i,j)   + match_cost;
            up_c   = C(i,j+1) + skip_cost;   
            left_c = C(i+1,j) + skip_cost;   

            [best, which] = min([diag_c, up_c, left_c]);
            C(i+1,j+1) = best;
            ptr(i+1,j+1) = which;
        end
    end

    i = n; j = m;
    pair_i = []; pair_j = [];
    while i > 0 || j > 0
        if i > 0 && j > 0 && ptr(i+1,j+1) == 1
            pair_i(end+1) = i; pair_j(end+1) = j; %#ok<AGROW>
            i = i - 1; j = j - 1;
        elseif i > 0 && (j == 0 || ptr(i+1,j+1) == 2)
            i = i - 1;
        else
            j = j - 1;
        end
    end
    pair_i = fliplr(pair_i); pair_j = fliplr(pair_j);

    used_e = false(n,1); used_i = false(m,1);
    for r = 1:numel(pair_i)
        ei = pair_i(r); ii = pair_j(r);
        used_e(ei) = true; used_i(ii) = true;
        matches(end+1) = struct('EmgIdx', emg_cand(ei).EmgIdx, 'ImuIdx', imu_cand(ii).ImuIdx, ...
            'EmgPre', emg_cand(ei).PreMarkerDur, 'EmgTrial', emg_cand(ei).TrialDur, 'EmgPost', emg_cand(ei).PostMarkerDur, ...
            'ImuPre', imu_cand(ii).PreMarkerDur, 'ImuTrial', imu_cand(ii).TrialDur, 'ImuPost', imu_cand(ii).PostMarkerDur, ...
            'EmgPreStartIdx', emg_cand(ei).PreStartIdx, 'EmgPostEndIdx', emg_cand(ei).PostEndIdx, ...
            'ImuPreStartIdx', imu_cand(ii).PreStartIdx, 'ImuPostEndIdx', imu_cand(ii).PostEndIdx, ...
            'EmgChronoIdx', emg_cand(ei).ChronoIdx, 'ImuChronoIdx', imu_cand(ii).ChronoIdx, 'ImuFileIndex', imu_cand(ii).FileIndex, 'InOrder', true); %#ok<AGROW>
    end

    unmatched_emg = emg_cand(~used_e);
    unmatched_imu = imu_cand(~used_i);
end


function [matches, unmatched_emg, unmatched_imu] = match_trials_unordered(emg_cand, imu_cand, tol)
    n = numel(emg_cand); m = numel(imu_cand);
    matches = struct('EmgIdx', {}, 'ImuIdx', {}, 'EmgPre', {}, 'EmgTrial', {}, 'EmgPost', {}, ...
        'ImuPre', {}, 'ImuTrial', {}, 'ImuPost', {}, 'EmgPreStartIdx', {}, 'EmgPostEndIdx', {}, ...
        'ImuPreStartIdx', {}, 'ImuPostEndIdx', {}, ...
        'EmgChronoIdx', {}, 'ImuChronoIdx', {}, 'ImuFileIndex', {}, 'InOrder', {});

    if n == 0 || m == 0
        unmatched_emg = emg_cand;
        unmatched_imu = imu_cand;
        return;
    end

    pairs = [];
    for i = 1:n
        for j = 1:m
            d_pre   = abs(emg_cand(i).PreMarkerDur - imu_cand(j).PreMarkerDur);
            d_trial = abs(emg_cand(i).TrialDur      - imu_cand(j).TrialDur);
            d_post  = abs(emg_cand(i).PostMarkerDur - imu_cand(j).PostMarkerDur);
            if d_pre <= tol && d_trial <= tol && d_post <= tol
                pairs = [pairs; i, j, d_pre + d_trial + d_post]; %#ok<AGROW>
            end
        end
    end

    used_e = false(n,1); used_i = false(m,1);
    if ~isempty(pairs)
        pairs = sortrows(pairs, 3);
        for r = 1:size(pairs,1)
            i = pairs(r,1); j = pairs(r,2);
            
            if used_i(j) || (used_e(i) && emg_cand(i).NumTrialsInRec <= 1)
                continue;
            end

            used_e(i) = true; used_i(j) = true;
            matches(end+1) = struct('EmgIdx', emg_cand(i).EmgIdx, 'ImuIdx', imu_cand(j).ImuIdx, ...
                'EmgPre', emg_cand(i).PreMarkerDur, 'EmgTrial', emg_cand(i).TrialDur, 'EmgPost', emg_cand(i).PostMarkerDur, ...
                'ImuPre', imu_cand(j).PreMarkerDur, 'ImuTrial', imu_cand(j).TrialDur, 'ImuPost', imu_cand(j).PostMarkerDur, ...
                'EmgPreStartIdx', emg_cand(i).PreStartIdx, 'EmgPostEndIdx', emg_cand(i).PostEndIdx, ...
                'ImuPreStartIdx', imu_cand(j).PreStartIdx, 'ImuPostEndIdx', imu_cand(j).PostEndIdx, ...
                'EmgChronoIdx', emg_cand(i).ChronoIdx, 'ImuChronoIdx', imu_cand(j).ChronoIdx, 'ImuFileIndex', imu_cand(j).FileIndex, 'InOrder', false); %#ok<AGROW>
        end
    end

    unmatched_emg = emg_cand(~used_e);
    unmatched_imu = imu_cand(~used_i);
end


function in_order = compute_lis_mask(seq)
    n = numel(seq);
    in_order = true(1, n);
    if n <= 1
        return;
    end
    dp = ones(1, n);
    parent = zeros(1, n);
    for i = 2:n
        for j = 1:i-1
            if seq(j) < seq(i) && dp(j) + 1 > dp(i)
                dp(i) = dp(j) + 1;
                parent(i) = j;
            end
        end
    end
    [~, idx] = max(dp);
    in_order = false(1, n);
    while idx > 0
        in_order(idx) = true;
        idx = parent(idx);
    end
end


function plot_sync_qc(t_vec, sig, segments, trials, label_str, out_png)
    max_pts = 20000;
    n = numel(t_vec);
    if n > max_pts
        step = ceil(n / max_pts);
        t_plot = t_vec(1:step:end);
        sig_plot = sig(1:step:end);
    else
        t_plot = t_vec; sig_plot = sig;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1400 400]);
    plot(t_plot, sig_plot, 'Color', [0.6 0.6 0.6]); hold on;
    yl = ylim;

    for tI = 1:numel(trials)
        tr = trials(tI);
        patch([tr.TrialStart, tr.TrialEnd, tr.TrialEnd, tr.TrialStart], ...
              [yl(1), yl(1), yl(2), yl(2)], [1.0 0.75 0.75], 'EdgeColor', 'none');
        text(mean([tr.TrialStart, tr.TrialEnd]), yl(2)*0.9, sprintf('trial %.2fs', tr.TrialDur), ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
        patch([tr.PreMarkerStart, tr.TrialStart, tr.TrialStart, tr.PreMarkerStart], ...
              [yl(1), yl(1), yl(2), yl(2)], [1.0 0.9 0.5], 'EdgeColor', 'none');
        text(mean([tr.PreMarkerStart, tr.TrialStart]), yl(2)*0.9, sprintf('%.2fs', tr.PreMarkerDur), ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
        patch([tr.TrialEnd, tr.PostMarkerEnd, tr.PostMarkerEnd, tr.TrialEnd], ...
              [yl(1), yl(1), yl(2), yl(2)], [1.0 0.9 0.5], 'EdgeColor', 'none');
        text(mean([tr.TrialEnd, tr.PostMarkerEnd]), yl(2)*0.9, sprintf('%.2fs', tr.PostMarkerDur), ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
    end

    plot(t_plot, sig_plot, 'Color', [0.3 0.3 0.3]);
    xlabel('Time (s)'); ylabel('Signal');
    title(label_str, 'Interpreter', 'none');
    hold off;
    print(fig, out_png, '-dpng', '-r100');
    close(fig);
end


function tbl = load_imu_table_v2(fpath)
    opts = detectImportOptions(fpath, 'FileType', 'text');
    tbl = readtable(fpath, opts);
end