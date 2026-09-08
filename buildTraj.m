%% AP
clear,clc 

axis_angles = [0.051523, -54.416751, 92.296563, 0.094545, -32.879845, -0.074910];
         
%-600
start = [     0.996195,    -0.000000,    -0.087156,  2462.560007 ;
             -0.000000,    -1.000000,     0.000000,    -1.827019 ;
             -0.087156,    -0.000000,    -0.996195,  1994.080988 ;
              0.000000,     0.000000,     0.000000,     1.000000 ];

radius = 1300;
roll = 0;
yaw = 0;

% 12deg/s vectors
%disp_names = ["disp_025_seg01.csv","disp_050_seg01.csv","disp_100.csv","disp_150.csv","disp_200.csv"];
%vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_100.csv", "angvel_150.csv", "angvel_200.csv"];

% 24deg/s vectors
% disp_names = "disp_24_100.csv";
% vel_names = "angvel_24_100.csv";

% Paired Vels
% 7deg/s
disp_names = ["disp_025_seg01.csv", "disp_050_seg01.csv", "disp_100.csv", "disp_150.csv", "disp_200.csv"];
vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_100.csv", "angvel_150.csv", "angvel_200.csv"];
% 0.5hz at 3.5deg/s (matching 0.25hz at 7deg/s
%disp_names = "disp_050_seg01.csv";
%vel_names = "angvel_050_seg01.csv";
% 2hz at 5.25deg/s
% disp_names = "disp_200.csv";
% vel_names = "angvel_200.csv";
% % 2hz at 7deg/s, but 15 points per cycle
%disp_names = "disp_150.csv";
%vel_names = "angvel_150.csv";
%disp_names = "disp_050_seg01.csv";
%vel_names = "angvel_050_seg01.csv";

% disp_names = ["disp_025_seg01.csv", "disp_050_seg01.csv","disp_200.csv"];
% vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_200.csv","angvel_200.csv"];

%axis = "sF";
axis = "FAx";

apo_dist = 5;


for i = 1:numel(disp_names)
    disp = readmatrix(disp_names(i));
    vel = readmatrix(vel_names(i));
    
    temp_name = string(disp_names(i));
    freq = regexp(temp_name, '_(\d+)', 'tokens', 'once');
    
    freq_val = str2double(freq) / 100; % Converts "025" string to 0.25
    fs = 1000;
    % Target a fixed resolution per cycle (e.g., 40 points per sine wave)
%     target_points_per_cycle = 15; 
%     raw_points_per_cycle = fs / freq_val;
%     down_sample = round(raw_points_per_cycle / target_points_per_cycle);
%     down_sample = max(1, down_sample); % Ensure it's at least 1
down_sample = 27; %number of samples per cycle
    
    filename = "AP_"+ freq + "_d"+down_sample+"_apo"+apo_dist+"_A50"+".csv";

    [position, orientation,~,~,tang_vel] = calcPose_Kuka(start, disp, roll, yaw, radius, vel);

%     pos_down1 = position(1:down_sample:end,:);
%     ori_down1 = orientation(1:down_sample:end,:);
    
    [pos_down, step] = downsample_matrix_and_dedup(position, down_sample, freq_val, 1000, 3);
    [ori_down,~] = downsample_matrix_and_dedup(orientation, down_sample, freq_val, 1000, 2);
    tang_vel_down1 = tang_vel(1:step:end,:);


% pos_down = downsample_multi_cycle(position, down_sample);
% ori_down = downsample_multi_cycle(orientation, down_sample);
% pos_down = phase_warped_downsample_multi(position,down_sample);
% ori_down = phase_warped_downsample_multi(orientation,down_sample);

     generateKUKALinearSRC(filename,pos_down,tang_vel_down1);
     generateKUKALinearDAT(filename,axis_angles,pos_down,ori_down,tang_vel_down1, apo_dist);
end


%% ML

axis_angles = [0.000340, -38.971005, 59.836584, -85.326283, 88.221407, 110.793006];

          
%-600
start = [    -0.000000,     0.999921,     0.012566,  2441.372992 ;
              0.996195,     0.001095,    -0.087149,   378.348007 ;
             -0.087156,     0.012518,    -0.996116,  2015.644016 ;
              0.000000,     0.000000,     0.000000,     1.000000 ];

          
down_sample = 25;
radius = 1300;
pitch = 0;
yaw = 0;

disp_names = ["disp_025_seg01.csv","disp_050_seg01.csv","disp_100.csv","disp_150.csv","disp_200.csv"];
vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_100.csv", "angvel_150.csv", "angvel_200.csv"];

for i = 1:numel(disp_names)
    disp = readmatrix(disp_names(i));
    vel = readmatrix(vel_names(i));
    
    temp_name = string(disp_names(i));
    freq = regexp(temp_name, '_(\d+)', 'tokens', 'once');
    
    freq_val = str2double(freq) / 100; % Converts "025" string to 0.25
    fs = 1000; %match fs used to generate the sine waves
    % Target a fixed resolution per cycle (e.g., 40 points per sine wave)
    target_points_per_cycle = 30; 
    raw_points_per_cycle = fs / freq_val;
    down_sample = round(raw_points_per_cycle / target_points_per_cycle);
    down_sample = max(1, down_sample); % Ensure it's at least 1
    
    filename = "ML_"+ freq + "_d"+down_sample+"_FAx.csv";

    [position, orientation,~,~,tang_vel] = calcPose_Kuka(start, pitch, disp, yaw, radius, vel);

%     pos_down1 = position(1:down_sample:end,:);
%     ori_down1 = orientation(1:down_sample:end,:);
    tang_vel_down1 = tang_vel(1:down_sample:end,:);
    pos_down = downsample_matrix_and_dedup(position, down_sample, 3);
    ori_down = downsample_matrix_and_dedup(orientation, down_sample, 3);

%     generateKUKALinearSRC(filename,pos_down,tang_vel_down1);
%     generateKUKALinearDAT(filename,axis_angles,pos_down,ori_down,tang_vel_down1);
end



%% Dia L

axis_angles = [0.108783, -42.773814, 72.780994, -67.816923, -48.436829, 53.603402];

%-600
start = [     0.714486,    -0.696852,    -0.062509,  2623.364196 ;
             -0.694200,    -0.717215,     0.060735,  -304.631080 ;
             -0.087156,    -0.000000,    -0.996195,  1932.196118 ;
              0.000000,     0.000000,     0.000000,     1.000000 ];

radius = 1300;
roll = 0;
yaw = 45;

disp_names = ["disp_025_seg01.csv","disp_050_seg01.csv","disp_100.csv","disp_150.csv","disp_200.csv"];
vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_100.csv", "angvel_150.csv", "angvel_200.csv"];

for i = 1:numel(disp_names)
    disp = readmatrix(disp_names(i));
    vel = readmatrix(vel_names(i));
    
    temp_name = string(disp_names(i));
    freq = regexp(temp_name, '_(\d+)', 'tokens', 'once');
    
    freq_val = str2double(freq) / 100; % Converts "025" string to 0.25
    fs = 1000;
    % Target a fixed resolution per cycle (e.g., 40 points per sine wave)
    target_points_per_cycle = 30; 
    raw_points_per_cycle = fs / freq_val;
    down_sample = round(raw_points_per_cycle / target_points_per_cycle);
    down_sample = max(1, down_sample); % Ensure it's at least 1
    
    filename = "DiaL_"+ freq + "_d"+down_sample+"_FAx.csv";

    [position, orientation,~,~,tang_vel] = calcPose_Kuka_diag_v3(start, disp, roll, yaw, radius, vel);

%     pos_down1 = position(1:down_sample:end,:);
%     ori_down1 = orientation(1:down_sample:end,:);
    tang_vel_down1 = tang_vel(1:down_sample:end,:);
    pos_down = downsample_matrix_and_dedup(position, down_sample, 3);
    ori_down = downsample_matrix_and_dedup(orientation, down_sample, 2);

    generateKUKALinearSRC(filename,pos_down,tang_vel_down1);
    generateKUKALinearDAT(filename,axis_angles,pos_down,ori_down,tang_vel_down1);
end


%%

%% Dia R

axis_angles = [-0.003720, -40.156924, 66.864986, -108.725064, 48.050672, 121.869055];

%-600
start = [    -0.000000,    -0.000000,     1.000000,  -600.000000 ;
              0.000000,    -1.000000,    -0.000000,     0.000000 ;
              1.000000,     0.000000,     0.000000,    90.000000 ;
              0.000000,     0.000000,     0.000000,     1.000000 ];

radius = 1300;
roll = 0;
yaw = -45;

disp_names = ["disp_025_seg01.csv","disp_050_seg01.csv","disp_100.csv","disp_150.csv","disp_200.csv"];
vel_names = ["angvel_025_seg01.csv", "angvel_050_seg01.csv", "angvel_100.csv", "angvel_150.csv", "angvel_200.csv"];

for i = 1:numel(disp_names)
    disp = readmatrix(disp_names(i));
    vel = readmatrix(vel_names(i));
    
    temp_name = string(disp_names(i));
    freq = regexp(temp_name, '_(\d+)', 'tokens', 'once');
    
    freq_val = str2double(freq) / 100; % Converts "025" string to 0.25
    fs = 1000;
    % Target a fixed resolution per cycle (e.g., 40 points per sine wave)
    target_points_per_cycle = 30; 
    raw_points_per_cycle = fs / freq_val;
    down_sample = round(raw_points_per_cycle / target_points_per_cycle);
    down_sample = max(1, down_sample); % Ensure it's at least 1
    
    filename = "DiaR_"+ freq + "_d"+down_sample+"_FAx.csv";

    [position, orientation,~,~,tang_vel] = calcPose_Kuka_diag_v3(start, disp, roll, yaw, radius, vel);

%     pos_down1 = position(1:down_sample:end,:);
%     ori_down1 = orientation(1:down_sample:end,:);
    tang_vel_down1 = tang_vel(1:down_sample:end,:);
    pos_down = downsample_matrix_and_dedup(position, down_sample, 3);
    ori_down = downsample_matrix_and_dedup(orientation, down_sample, 2);

    generateKUKALinearSRC(filename,pos_down,tang_vel_down1);
    generateKUKALinearDAT(filename,axis_angles,pos_down,ori_down,tang_vel_down1);
end



%% --------------------------------------

%----------------------------------------
%% Main function: Downsample rows + handle rounded duplicates on specified column
% function data_clean = downsample_matrix_and_dedup(data, decimation_factor, col_to_check)
%     % data: n x 3 matrix
%     % decimation_factor: e.g. 40
%     % col_to_check: integer (1, 2, or 3) - column to check for duplicates
%     
%     if nargin < 3
%         col_to_check = 1;  % default to first column
%     end
%     
%     if col_to_check < 1 || col_to_check > size(data,2)
%         error('col_to_check must be a valid column index (1 to %d)', size(data,2));
%     end
%     
%     % Step 1: Downsample rows (keep every decimation_factor-th row)
%     data_down = data(1:decimation_factor:end, :);
%     
%     % Step 2: Clean duplicates only on the specified column
%     data_clean = remove_rounded_duplicates_specific_col(data_down, col_to_check);
%     
%     fprintf('Original rows: %d\n', size(data,1));
%     fprintf('Downsampled rows: %d\n', size(data_clean,1));
%     fprintf('Duplicate checking applied to column: %d\n', col_to_check);
% end

function [data_clean,step] = downsample_matrix_and_dedup(data, samples_per_cycle, freq_hz, fs, col_to_check)
    % data: n x 3 (or more) matrix
    % samples_per_cycle: desired number of points per cycle (e.g. 8, 10, 12)
    % freq_hz: exact frequency of the original sine wave in Hz (e.g. 0.25 or 2.0)
    % fs: original sampling frequency in Hz (e.g. 4400, 10000, etc.)
    % col_to_check: column to check for consecutive duplicates (default = 1)
    
    if nargin < 5 || isempty(col_to_check)
        col_to_check = 1;
    end
    if nargin < 4 || isempty(fs)
        error('You must provide the original sampling frequency fs');
    end
    
    L = size(data,1);
    
    % Calculate step size
    samples_per_sec_desired = samples_per_cycle * freq_hz;
    step = max(1, round(fs / samples_per_sec_desired));
    
    % Downsample
    data_down = data(1:step:end, :);
    
    % Remove consecutive duplicates
    data_clean = remove_rounded_duplicates_specific_col(data_down, col_to_check);
    
    % Info
    fprintf('Original rows: %d → Downsampled rows: %d\n', L, size(data_clean,1));
    fprintf('Samples per cycle: %d | Signal freq: %.3f Hz | Original fs: %.1f Hz | step: %.3f\n', ...
            samples_per_cycle, freq_hz, fs, step);
end

%% Core deduplication - only on one specified column
function y = remove_rounded_duplicates_specific_col(x, col)
    y = x;                    % Work on copy
    n_rows = size(x, 1);
    
    fprintf('Processing column %d for duplicates...\n', col);
    i = 2;
    while i < n_rows
        % Round to 2 decimal places for duplicate check (only on chosen column)
        rnd_prev = round(y(i-1, col), 2);
        rnd_curr = round(y(i,   col), 2);
        
        if abs(rnd_prev - rnd_curr) < 1e-9
            % Duplicate found -> average prev and next in THIS column only
            if i < n_rows
                avg_value = (y(i-1, col) + y(i+1, col)) / 2;
                y(i, col) = avg_value;
                fprintf('  Fixed duplicate in col %d at row %d: %.4f → %.4f\n', ...
                        col, i, x(i,col), y(i,col));
            else
                % Last row
                y(i, col) = (y(i-1, col) + y(i, col)) / 2;
            end
            i = i + 1;
        else
            i = i + 1;
        end
    end
end


%%
function downsampled = downsample_multi_cycle(data, target_N)
    % data: NxM matrix
    % target_N: desired number of rows after downsampling
    
    L = size(data, 1);
    if L <= target_N
        downsampled = data;
        return;
    end
    
    signal = data(:,1);                    % use first column as reference
    
    % Find peaks
    signal_range = max(signal) - min(signal);
    [~, peak_locs] = findpeaks(signal, 'MinPeakProminence', 0.4*signal_range);
    
    num_cycles_approx = max(1, length(peak_locs));
    
    % Warping
    t = linspace(0, 2*pi*num_cycles_approx, target_N);
    warped = abs(sin(t));                  
    warped = warped / max(warped);
    
    cum_warped = cumsum(warped) / sum(warped);
    
    idx = round(1 + (L-1) * cum_warped);
    idx = unique(max(1, min(L, idx)));
    
    if length(idx) > target_N
        idx = idx(round(linspace(1, length(idx), target_N)));
    end
    
    downsampled = data(idx, :);
end

%%

function downsampled = phase_warped_downsample_multi(data, target_N)
    % data: NxM matrix
    % target_N: TOTAL number of points you want after downsampling
    
    L = size(data, 1);
    if L <= target_N
        downsampled = data;
        return;
    end
    
    % Reference signal (first column)
    signal = data(:,1) - mean(data(:,1));
    
    % Detect number of cycles
    signal_range = max(signal) - min(signal);
    [~, peak_locs] = findpeaks(signal, 'MinPeakProminence', 0.3*signal_range);
    [~, trough_locs] = findpeaks(-signal, 'MinPeakProminence', 0.3*signal_range);
    
    num_cycles = max(1, max(length(peak_locs), length(trough_locs)));
    
    % Warping across all cycles
    t = linspace(0, 2*pi * num_cycles, target_N);
    warped = abs(sin(t));                  
    warped = warped / max(warped);
    
    cum_warped = cumsum(warped) / sum(warped);
    
    idx = round(1 + (L-1) * cum_warped);
    idx = unique(max(1, min(L, idx)));
    
    if length(idx) > target_N
        idx = idx(round(linspace(1, length(idx), target_N)));
    end
    
    downsampled = data(idx, :);
end