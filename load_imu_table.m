function data = load_imu_table(filename)
    % LOAD_IMU_TABLE Loads IMU data from text/CSV files and fixes the missing Sync column
    
    % Try to find the file on the MATLAB path
    fullpath = which(filename);
    if isempty(fullpath)
        if isfile(filename)
            fullpath = fullfile(pwd, filename);
        else
            error('File not found: %s\nCheck that the file exists and the folder is on the MATLAB path.', filename);
        end
    end
    
    % Create import options
    opts = detectImportOptions(fullpath);
    opts.Delimiter = ',';
    opts.VariableNamesLine = 1;
    opts.DataLines = [2 Inf];
    
    % Read the table
    data = readtable(fullpath, opts);
    data.Properties.VariableNames = strtrim(data.Properties.VariableNames);
    
    w = width(data);
    
    %% Fix 'Sync' designation and align IMU 1 headers dynamically
    if w >= 10
        current_names = data.Properties.VariableNames;
        
        % 1. Retain the first 9 headers: Time (1), AI0-AI7 (2:9)
        corrected_names = current_names(1:9);
        
        % 2. Force the 10th column to be our synchronization axis
        corrected_names{10} = 'Sync';
        
        % 3. Align IMU 1 variables dynamically based on how many columns exist
        if w >= 11, corrected_names{11} = 'Ax1'; end
        if w >= 12, corrected_names{12} = 'Ay1'; end
        if w >= 13, corrected_names{13} = 'Az1'; end
        if w >= 14, corrected_names{14} = 'Gx1'; end
        if w >= 15, corrected_names{15} = 'Gy1'; end
        if w >= 16, corrected_names{16} = 'Gz1'; end
        
        % 4. Append whatever trailing columns remain (like IMU 2 arrays) unchanged
        if w > 16
            corrected_names(17:w) = current_names(17:end);
        end
        
        data.Properties.VariableNames = corrected_names;
    end
    
    fprintf('Loaded %d rows x %d columns from %s\n', ...
            height(data), width(data), filename);
end