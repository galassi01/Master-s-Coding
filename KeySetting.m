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
    calib1_file = 'KUKA_neckEMG_P04_AllOri1.txt';
    calib2_file = 'KUKA_neckEMG_P04_AllOri2.txt';
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

outfile = sprintf('P%02d_sync_summary', participant_num);
sync_out_dir = fullfile(pwd, outfile);   % where the matcher script wrote its CSVs
matched_csv  = fullfile(sync_out_dir, 'EMG_IMU_Matched_Trials.csv');
key_csv      = fullfile(sync_out_dir, 'IMU_Key_Table.csv');
keyed_csv    = fullfile(sync_out_dir, 'EMG_IMU_Matched_Trials_Keyed.csv');

if ~exist(sync_out_dir, 'dir')
    mkdir(sync_out_dir);
    fprintf('Created output folder: %s\n', sync_out_dir);
end

valid_conditions = {'EC', 'EO'};

% --- Reference lookup: (ID, Group label) -> code, used by the bulk block-
%     table entry mode (mode [2]). e.g. ID=7, Group="AP" -> "0800_0360".
%     Edit this table directly if the reference values change. ---
code_lookup_data = {
    1, 'ML',   '1320_0402';
    2, 'ML',   '0800_0393';
    3, 'ML',   '0550_0400';
    4, 'ML',   '0410_0408';
    5, 'ML',   '1320_0300';
    6, 'ML',   '0800_0200';
    7, 'ML',   '0550_0200';
    1, 'DiaR', '1330_0413';
    2, 'DiaR', '0830_0401';
    3, 'DiaR', '0480_0403';
    4, 'DiaR', '0370_0405';
    5, 'DiaR', '0390_0300';
    1, 'DiaL', '1330_0413';
    2, 'DiaL', '0830_0401';
    3, 'DiaL', '0480_0403';
    4, 'DiaL', '0370_0405';
    5, 'DiaL', '0390_0300';
    1, 'AP',   '1230_0695';
    2, 'AP',   '0990_0679';
    3, 'AP',   '0840_0679';
    4, 'AP',   '0380_0679';
    5, 'AP',   '1250_0290';
    6, 'AP',   '1000_0360';
    7, 'AP',   '0800_0360';
    8, 'AP',   '0410_0360';
};
code_lookup = cell2table(code_lookup_data, 'VariableNames', {'ID', 'Group', 'Code'});
code_lookup.Group = string(code_lookup.Group);
code_lookup.Code  = string(code_lookup.Code);

%% 1. Load or initialize the key table ====================================
imu_files  = dir(fullfile(imu_folder, 'KUKA_*.txt'));
file_names = string({imu_files.name}');

if isempty(file_names) && ~(exist(key_csv, 'file'))
    error(['dir(fullfile(''%s'', ''KUKA_*.txt'')) found ZERO files, and no existing key_csv to fall back on.\n' ...
        'This is almost certainly why every bulk-tag lookup fails - the key table would be empty.\n' ...
        'Check that imu_folder is correct and that this machine can actually see that path right now\n' ...
        '(e.g. a mounted drive not currently mounted in this session would cause exactly this).'], imu_folder);
elseif isempty(file_names)
    fprintf(['WARNING: dir(fullfile(''%s'', ''KUKA_*.txt'')) found ZERO files. Falling back to whatever is\n' ...
        'already in %s, but if that''s stale/incomplete, bulk-tag lookups for real files will fail.\n' ...
        'Double check imu_folder is reachable before relying on this run.\n'], imu_folder, key_csv);
end

if exist(key_csv, 'file')
    key_table = readtable(key_csv, 'TextType', 'string');
    % Backward-compat: older key tables saved before the lookup-based bulk
    % mode existed won't have these columns yet - add them blank, BEFORE
    % any concatenation below (table vertcat needs matching columns).
    if ~ismember('GroupID', key_table.Properties.VariableNames)
        key_table.GroupID = nan(height(key_table), 1);
    end
    if ~ismember('GroupLabel', key_table.Properties.VariableNames)
        key_table.GroupLabel = repmat("", height(key_table), 1);
    end

    % Add any files that exist on disk but aren't in the key table yet
    % (e.g. new files added since the last time this was run)
    missing = setdiff(file_names, key_table.FileName);
    if ~isempty(missing)
        add_tbl = table(missing, repmat("", numel(missing), 1), repmat("", numel(missing), 1), ...
            nan(numel(missing), 1), repmat("", numel(missing), 1), ...
            'VariableNames', {'FileName', 'ParameterSet', 'Condition', 'GroupID', 'GroupLabel'});
        key_table = [key_table; add_tbl]; %#ok<AGROW>
        fprintf('Added %d new file(s) found on disk to the key table (unassigned).\n', numel(missing));
    end
else
    key_table = table(file_names, repmat("", numel(file_names), 1), repmat("", numel(file_names), 1), ...
        nan(numel(file_names), 1), repmat("", numel(file_names), 1), ...
        'VariableNames', {'FileName', 'ParameterSet', 'Condition', 'GroupID', 'GroupLabel'});
    fprintf('No existing key table found - starting a new one with %d file(s).\n', numel(file_names));
end

% Sort into natural (block, trial) NUMERIC order rather than whatever
% order dir() / the CSV happened to be in. dir() sorts filenames as text,
% so "block1_10" sorts right after "block1_1" and before "block1_2" -
% block1_1, block1_10, block1_11, ..., block1_19, block1_2, block1_20, ...
% That's confusing (and error-prone) when you're picking rows by number
% in mode [1], so the key table is always kept in true 1, 2, 3, ..., 9,
% 10, 11, ... order per block instead.
block_nums = nan(height(key_table), 1);
trial_nums = nan(height(key_table), 1);
for fi = 1:height(key_table)
    [block_nums(fi), trial_nums(fi)] = parse_block_trial(key_table.FileName(fi));
end
[~, sort_ix] = sortrows([block_nums, trial_nums], [1 2]);  % NaNs (unparseable names) sort last
key_table = key_table(sort_ix, :);

%% 2. Interactive tagging loop =============================================
fprintf('\n=== IMU trial key entry ===\n');

while true
    n_unassigned = sum(key_table.ParameterSet == "");
    fprintf('\n%d of %d files still unassigned.\n', n_unassigned, height(key_table));
    fprintf('Choose an entry mode:\n');
    fprintf('  [1] Tag a file selection with ONE ParameterSet + Condition (patterns or row numbers)\n');
    fprintf('  [2] Paste a block table (one row per trial, in trial order - ID + Group looked up to a code)\n');
    fprintf('  [q] Done\n');
    mode = strtrim(lower(input('Mode: ', 's')));

    switch mode
        case '1'
            param = strtrim(string(input('Parameter set (e.g. 05_4): ', 's')));
            if param == ""
                continue;
            end
            cond = "";
            while ~any(strcmpi(cond, valid_conditions))
                cond = strtrim(string(input('Condition (EC or EO): ', 's')));
                cond = upper(cond);
                if ~any(strcmpi(cond, valid_conditions))
                    fprintf('  Condition must be EC or EO - try again.\n');
                end
            end

            fprintf('\nCurrent key table:\n');
            disp(key_table);
            fprintf(['\nWhich files get ParameterSet="%s", Condition="%s"?\n' ...
                '  - row numbers, comma-separated (e.g. "3,4,5")\n' ...
                '  - filename substring(s), comma-separated, quotes optional\n' ...
                '      e.g.  block6_1, block6_6\n' ...
                '      or    "block6_1", "block6_6"\n' ...
                '  - "all" for every currently-unassigned file\n'], param, cond);
            sel = input('Selection: ', 's');

            idx = resolve_selection(sel, key_table);
            if isempty(idx)
                fprintf('  No files matched that selection - nothing changed.\n');
                continue;
            end

            key_table.ParameterSet(idx) = param;
            key_table.Condition(idx)    = cond;
            writetable(key_table, key_csv);
            fprintf('  Tagged %d file(s):\n', numel(idx));
            disp(key_table(idx, :));

        case '2'
            key_table = bulk_tag_block_table(key_table, code_lookup);
            writetable(key_table, key_csv);

        case 'q'
            break;

        otherwise
            fprintf('  Unrecognized option - enter 1, 2, or q.\n');
    end
end

writetable(key_table, key_csv);
fprintf('\nFinal key table saved to: %s\n', key_csv);
n_unassigned = sum(key_table.ParameterSet == "");
if n_unassigned > 0
    fprintf('NOTE: %d file(s) still have no ParameterSet/Condition assigned:\n', n_unassigned);
    disp(key_table(key_table.ParameterSet == "", :));
end

%% 3. Merge the key table into the matched-trials output ==================
if ~exist(matched_csv, 'file')
    fprintf('\nNOTE: %s not found - run the sync extractor/matcher script first if you want the merged output.\n', matched_csv);
    fprintf('The key table above has been saved regardless, and can be merged in later.\n');
    return;
end

match_table = readtable(matched_csv, 'TextType', 'string');
if ~ismember('IMU_FileName', match_table.Properties.VariableNames)
    error('Expected a column named "IMU_FileName" in %s - check the matcher script output format.', matched_csv);
end

merged = outerjoin(match_table, key_table, ...
    'LeftKeys', 'IMU_FileName', 'RightKeys', 'FileName', ...
    'MergeKeys', true, 'Type', 'left');

writetable(merged, keyed_csv);
fprintf('\nMerged keyed + matched table saved to: %s\n', keyed_csv);

n_unkeyed_matches = sum(merged.ParameterSet == "" | ismissing(merged.ParameterSet));
if n_unkeyed_matches > 0
    fprintf('NOTE: %d matched trial row(s) have no ParameterSet/Condition (their IMU file wasn''t tagged yet).\n', n_unkeyed_matches);
end

%% 4. Example filters ======================================================
fprintf('\n=== Example filters on the keyed table ===\n');
fprintf('%% All EC trials from parameter set 05_4:\n');
fprintf('  filt = merged(merged.Condition == "EC" & merged.ParameterSet == "05_4", :);\n');
fprintf('%% All trials for a given parameter set, either condition:\n');
fprintf('  filt = merged(merged.ParameterSet == "05_4", :);\n');
fprintf('%% All EO trials across every parameter set:\n');
fprintf('  filt = merged(merged.Condition == "EO", :);\n');
fprintf('%% Count trials per (ParameterSet, Condition) combination:\n');
fprintf('  summary_counts = groupsummary(merged, {''ParameterSet'',''Condition''});\n');

if exist('merged', 'var')
    example_filt = merged(merged.Condition == "EC", :); %#ok<NASGU>
    fprintf('\n(You now have `merged` in the workspace to filter directly, e.g. try the lines above.)\n');
end


%% ========================================================================
%  HELPER FUNCTIONS
% ========================================================================

function idx = resolve_selection(sel, key_table)
% RESOLVE_SELECTION  Turn a user's typed selection into row indices into
% key_table. Supports:
%   - "all"                              -> every currently-unassigned row
%   - "3,4,5" / "3, 7, 12"               -> those row numbers directly
%   - one or more filename substrings, comma-separated, quotes optional,
%     matched as an OR (a file matching ANY of them is included):
%       block6_1, block6_6
%       "block6_1", "block6_6"
% This was the bug in the previous version: a multi-pattern selection like
% "block6_1", "block6_6" (including the literal quotes and comma) was
% passed to `contains` as ONE string, which no filename literally contains,
% so nothing ever matched. Now each comma-separated piece (with any
% surrounding quotes stripped) is matched independently.

    sel = strtrim(sel);
    n = height(key_table);

    if strcmpi(sel, 'all')
        idx = find(key_table.ParameterSet == "");
        return;
    end

    % Try parsing as a comma-separated list of row numbers first
    parts = strsplit(sel, ',');
    parts_trimmed = strtrim(parts);
    nums = str2double(parts_trimmed);
    if all(~isnan(nums)) && ~isempty(nums)
        nums = round(nums);
        valid = nums >= 1 & nums <= n;
        if any(~valid)
            fprintf('  Some row numbers were out of range (1-%d) and were ignored.\n', n);
        end
        idx = nums(valid);
        idx = idx(:);
        return;
    end

    % Otherwise: one or more filename substrings, comma-separated, with
    % optional surrounding double or single quotes on each piece
    patterns = cell(1, numel(parts_trimmed));
    for p = 1:numel(parts_trimmed)
        piece = parts_trimmed{p};
        piece = regexprep(piece, '^["'']|["'']$', '');  % strip leading/trailing quote chars
        patterns{p} = strtrim(piece);
    end
    patterns = patterns(~cellfun(@isempty, patterns));

    mask = false(n, 1);
    for p = 1:numel(patterns)
        mask = mask | contains(lower(key_table.FileName), lower(patterns{p}));
    end
    idx = find(mask);
end


function key_table = bulk_tag_block_table(key_table, code_lookup)
% BULK_TAG_BLOCK_TABLE  Paste a small table for one block, one row per
% trial IN TRIAL ORDER: row 1 -> trial 1 of that block, row 2 -> trial 2,
% etc. Each line has exactly 3 tokens: <ID> <Group> <Condition>, e.g.
%
%   7  AP  EO
%   4  AP  EC
%   3  AP  EO
%
% The Condition (last token) must be EC/EO. The ID + Group (first two
% tokens) are NOT concatenated into the ParameterSet directly - they're a
% lookup key into `code_lookup` (ID, Group) -> Code, e.g. ID=7, Group="AP"
% looks up to "0800_0360", and THAT code becomes the ParameterSet. The raw
% ID/Group are also stored in the GroupID/GroupLabel columns for
% traceability, in case you need to check which row produced which code.
%
% This matches to files via the "block<N>_<trial>" pattern already in the
% filename (parse_block_trial), so it only tags files that actually exist
% for that block/trial combination.

    valid_conditions = {'EC', 'EO'};

    block_str = input('Block number (e.g. 4 for "block4_..."): ', 's');
    block_num = str2double(block_str);
    if isnan(block_num)
        fprintf('  Not a valid block number - cancelled.\n');
        return;
    end

    % --- Diagnostic: show what's actually available for this block BEFORE
    %     attempting to match individual rows. If this comes back empty,
    %     the problem is upstream (key_table doesn't have this block's
    %     files at all - wrong imu_folder, stale key_csv, etc.) rather
    %     than anything about the pasted rows themselves. ---
    all_blk = nan(height(key_table), 1);
    all_trl = nan(height(key_table), 1);
    for fi = 1:height(key_table)
        [all_blk(fi), all_trl(fi)] = parse_block_trial(key_table.FileName(fi));
    end
    this_block_mask = all_blk == block_num;
    fprintf('Key table has %d file(s) total; %d of them parsed as block %d.\n', ...
        height(key_table), sum(this_block_mask), block_num);
    if any(this_block_mask)
        fprintf('Their trial numbers are: %s\n', mat2str(sort(all_trl(this_block_mask))'));
    else
        fprintf(['WARNING: no files in the key table parsed as block %d at all - every row below will fail\n' ...
            'to find a match. This means either the key table is missing this block''s files entirely\n' ...
            '(check imu_folder / re-run Section 1), or the block number you entered is wrong.\n' ...
            'A few filenames currently in the key table, for reference:\n'], block_num);
        sample_n = min(5, height(key_table));
        if sample_n == 0
            fprintf('  (key table is completely empty)\n');
        else
            for si = 1:sample_n
                fprintf('  %s\n', key_table.FileName(si));
            end
        end
    end

    fprintf(['Paste rows now, one trial per line, IN TRIAL ORDER (row 1 = trial 1, row 2 = trial 2, ...).\n' ...
        'Each line: <ID> <Group> <Condition>   e.g.  7 AP EO   (spaces, commas, tabs, or\n' ...
        'semicolons all work as separators - "7, AP, EO" and "7 AP EO" are both fine)\n' ...
        'Enter a blank line when done.\n']);

    rows = {};
    while true
        ln = input('', 's');
        if isempty(strtrim(ln))
            break;
        end
        rows{end+1} = ln; %#ok<AGROW>
    end

    n_rows = numel(rows);
    if n_rows == 0
        fprintf('  No rows entered - nothing changed.\n');
        return;
    end

    n_tagged = 0;
    for trial_num = 1:n_rows
        tokens = strsplit(strtrim(rows{trial_num}), {' ', '\t', ',', ';'});
        tokens = tokens(~cellfun(@isempty, tokens));
        if numel(tokens) ~= 3
            fprintf('  Row %d ("%s") needs exactly 3 tokens (ID, Group, Condition) - got %d, skipped.\n', ...
                trial_num, rows{trial_num}, numel(tokens));
            continue;
        end

        id_val   = str2double(tokens{1});
        group_lbl = tokens{2};
        cond     = upper(tokens{3});

        if isnan(id_val)
            fprintf('  Row %d: "%s" is not a valid numeric ID - skipped.\n', trial_num, tokens{1});
            continue;
        end
        if ~any(strcmpi(cond, valid_conditions))
            fprintf('  Row %d: last token "%s" is not EC/EO - skipped.\n', trial_num, tokens{3});
            continue;
        end

        lookup_row = code_lookup.ID == id_val & strcmpi(code_lookup.Group, group_lbl);
        if ~any(lookup_row)
            fprintf('  Row %d: no lookup entry for ID=%g, Group="%s" - skipped.\n', trial_num, id_val, group_lbl);
            continue;
        end
        code = code_lookup.Code(find(lookup_row, 1));

        target_idx = find(all_blk == block_num & all_trl == trial_num, 1);

        if isempty(target_idx)
            fprintf('  Row %d: no file found for block%d_%d - skipped.\n', trial_num, block_num, trial_num);
            continue;
        end

        key_table.ParameterSet(target_idx) = code;
        key_table.Condition(target_idx)    = string(cond);
        key_table.GroupID(target_idx)      = id_val;
        key_table.GroupLabel(target_idx)   = string(group_lbl);
        n_tagged = n_tagged + 1;
        fprintf('  block%d_%d -> ParameterSet="%s" (from ID=%g, Group="%s"), Condition="%s"\n', ...
            block_num, trial_num, code, id_val, group_lbl, cond);
    end

    fprintf('Tagged %d of %d pasted row(s).\n', n_tagged, n_rows);
end


function [blockNum, trialNum] = parse_block_trial(fname)
% PARSE_BLOCK_TRIAL  Extract the block and trial number from a filename
% like "KUKA_neckEMG_P04_Trial_Block1_1.txt" -> blockNum=1, trialNum=1.
% Case-insensitive on purpose - actual filenames use "Block" (capital B),
% and a case-sensitive regexp here silently matched nothing at all,
% which is why every bulk-tag lookup was failing uniformly.
    fname = char(fname);
    tok = regexpi(fname, 'block(\d+)_(\d+)', 'tokens', 'once');
    if isempty(tok)
        blockNum = NaN; trialNum = NaN;
    else
        blockNum = str2double(tok{1});
        trialNum = str2double(tok{2});
    end
end
