% --- Configuration ---
Fs = 2000;                  % Sampling rate (Hz)
windowDuration = 0.10;     % 100 ms window
windowLen = round(windowDuration * Fs); % 200 samples
est_freq = 1.32;

% --- Load / Input Rectified Signal ---
% Assume 'rEMG' is your input rectified vector
EMG1 = data1.emg.signal_highpass_cycletrim;
EMG2 = data2.emg.signal_highpass_cycletrim;

IMU1 = data1.imu_flange.six_axis_resamp;
IMU2 = data2.imu_flange.six_axis_resamp;


% --- 1. Centered 100ms Moving RMS (Recommended for Offline Analysis) ---
% Preserves timing alignment with movement kinematics/force without phase delay.
emg_rms1 = sqrt(movmean(EMG1.^2, windowLen,2));
emg_rms2 = sqrt(movmean(EMG2.^2, windowLen,2));

[channel_means, rms_cycles, cycle_lengths, outlier_cycles, raw_cycle_lengths, median_len] = segment_and_average_emg(emg_rms1, IMU1, est_freq, 'direction_axis',4);
fprintf('Expected: %.1f samples | Actual range: [%d, %d] | Median: %.1f\n', ...
    expected_period_samples, min(raw_cycle_lengths), max(raw_cycle_lengths), median(raw_cycle_lengths));
% [channel_means, rms_cycles, cycle_lengths, median_lenm crossings_emg]= segment_and_average_emg(emg_rms1,IMU1, 'direction_axis',4);
% 

% function [channel_mean_rms, rms_cycles_cell, raw_cycle_lengths, median_len, crossings_emg] = segment_and_average_emg(EMG_rms, imu_signal, varargin)
% % SEGMENT_AND_AVERAGE_EMG Segments a 64xN continuous RMS EMG matrix into cycles.
% % Interpolates each natural cycle to the MEDIAN cycle length (L_median) across 
% % all 64 channels so that every cell has identical dimensions [1 x L_median].
% %
% % INPUTS:
% %   EMG_rms    : [64 x N_emg] continuous RMS matrix (e.g., EMG_rms1)
% %   imu_signal : [N_imu x 6] or [N_imu x 1] IMU matrix/vector
% %
% % OUTPUTS:
% %   channel_mean_rms  : [64 x L_median] Cycle-averaged RMS matrix per channel
% %   rms_cycles_cell   : {64 x N_cycles} Cell array where EVERY cell is [1 x L_median]
% %   raw_cycle_lengths : [1 x N_cycles] Natural, un-interpolated sample length per cycle
% %   median_len        : Target median cycle length (in samples)
% %   crossings_emg     : Detected cycle boundary indices in EMG domain
% 
% %% --- Parse Optional Inputs ---
% p = inputParser;
% addParameter(p, 'direction_axis', 4, @(x) isnumeric(x) && x >= 1);
% addParameter(p, 'cycle_type', 'full', @(x) ischar(x) || isstring(x));
% parse(p, varargin{:});
% 
% rot_col    = p.Results.direction_axis;
% cycle_type = lower(char(p.Results.cycle_type));
% 
% %% --- Orient Input Matrices ---
% if size(EMG_rms, 1) ~= 64 && size(EMG_rms, 2) == 64
%     EMG_rms = EMG_rms';
% end
% if size(imu_signal, 1) < size(imu_signal, 2)
%     imu_signal = imu_signal';
% end
% if size(imu_signal, 2) >= rot_col
%     imu_vec = imu_signal(:, rot_col);
% else
%     imu_vec = imu_signal(:, 1);
% end
% 
% %% --- 1. Detect Debounced Zero-Crossings ---
% imu_centered = imu_vec - mean(imu_vec);
% raw_zc = find(imu_centered(1:end-1) .* imu_centered(2:end) <= 0);
% 
% if isempty(raw_zc)
%     error('No zero crossings detected in the provided IMU signal.');
% end
% 
% median_half_period = median(diff(raw_zc));
% debounce_samples   = max(1, round(0.25 * median_half_period));
% 
% crossings_imu = raw_zc(1);
% for k = 2:length(raw_zc)
%     if (raw_zc(k) - crossings_imu(end)) >= debounce_samples
%         crossings_imu = [crossings_imu; raw_zc(k)]; %#ok<AGROW>
%     end
% end
% 
% %% --- 2. Map Crossings to EMG & Form Cycle Boundaries ---
% N_imu = length(imu_vec);
% N_emg = size(EMG_rms, 2);
% crossings_emg = 1 + round((crossings_imu - 1) * ((N_emg - 1) / (N_imu - 1)));
% 
% if strcmp(cycle_type, 'full')
%     valid_starts = 1:2:(length(crossings_emg) - 2);
%     cycle_starts = crossings_emg(valid_starts);
%     cycle_ends   = crossings_emg(valid_starts + 2);
% else % 'half'
%     cycle_starts = crossings_emg(1 : end-1);
%     cycle_ends   = crossings_emg(2 : end);
% end
% 
% num_cycles = length(cycle_starts);
% if num_cycles < 1
%     error('Not enough zero-crossings to form at least one complete cycle.');
% end
% 
% %% --- 3. Compute Natural Lengths & Target Median Length ---
% raw_cycle_lengths = (cycle_ends - cycle_starts) + 1;
% median_len        = round(median(raw_cycle_lengths));
% 
% fprintf('\n--- MEDIAN CYCLE INTERPOLATION REPORT ---\n');
% fprintf('Total Cycles Detected : %d\n', num_cycles);
% fprintf('Target Median Length  : %d samples\n', median_len);
% fprintf('Min Natural Length    : %d samples\n', min(raw_cycle_lengths));
% fprintf('Max Natural Length    : %d samples\n', max(raw_cycle_lengths));
% fprintf('Mean ± Std            : %.2f ± %.2f samples\n', mean(raw_cycle_lengths), std(raw_cycle_lengths));
% fprintf('-----------------------------------------\n\n');
% 
% %% --- 4. Interpolate Each Cycle to Median Length ---
% rms_cycles_cell = cell(64, num_cycles);
% target_grid     = linspace(0, 1, median_len);
% 
% for c = 1:num_cycles
%     i_start = cycle_starts(c);
%     i_end   = cycle_ends(c);
% 
%     % Slice raw cycle data across all 64 channels
%     raw_slice     = EMG_rms(:, i_start:i_end); % [64 x N_raw_samples]
%     n_raw_samples = size(raw_slice, 2);
%     orig_grid     = linspace(0, 1, n_raw_samples);
% 
%     % Interpolate each channel to median_len using pchip
%     for ch = 1:64
%         rms_cycles_cell{ch, c} = interp1(orig_grid, raw_slice(ch, :), target_grid, 'pchip');
%     end
% end
% 
% %% --- 5. Compute Channel Means ---
% channel_mean_rms = zeros(64, median_len);
% for ch = 1:64
%     % Stack row-wise: [N_cycles x L_median]
%     ch_cycles_matrix = vertcat(rms_cycles_cell{ch, :});
%     channel_mean_rms(ch, :) = mean(ch_cycles_matrix, 1);
% end
% 
% end


%%

function [channel_mean_rms, rms_cycles_cell, num_cycles, outlier_cycles, raw_cycle_lengths, median_len] = segment_and_average_emg(EMG_rms, imu_signal, input_freq_hz, varargin)
% SEGMENT_AND_AVERAGE_EMG Segments a 64xN continuous RMS EMG matrix into cycles.
% Flags and excludes outlier cycles deviating >5% from the expected cycle period,
% then interpolates valid cycles to MEDIAN cycle length (L_median) across 64 channels.
%
% INPUTS:
%   EMG_rms       : [64 x N_emg] continuous RMS matrix (e.g., EMG_rms1)
%   imu_signal    : [N_imu x 6] or [N_imu x 1] IMU matrix/vector
%   input_freq_hz : Expected fundamental movement frequency in Hz (e.g., 1.32 Hz)
%
% OPTIONAL PARAMETERS (Name-Value):
%   'emg_fs'         : EMG sampling rate in Hz (Default: 2000)
%   'direction_axis' : Column of IMU to use for segmentation (Default: 4 for Gx)
%   'cycle_type'     : 'full' (3 crossings per cycle) or 'half' (2 crossings) (Default: 'full')
%
% OUTPUTS:
%   channel_mean_rms  : [64 x L_median] Cycle-averaged RMS matrix per channel (valid cycles only)
%   rms_cycles_cell   : {64 x N_valid_cycles} Cell array of interpolated valid cycles [1 x L_median]
%   num_cycles        : Total number of initial cycles detected
%   outlier_cycles    : Vector of cycle indices flagged and excluded due to >5% period deviation
%   raw_cycle_lengths : [1 x N_cycles] Natural sample length per cycle (all initial cycles)
%   median_len        : Target median cycle length (in samples across valid cycles)

%% --- Parse Inputs ---
p = inputParser;
addRequired(p, 'EMG_rms', @isnumeric);
addRequired(p, 'imu_signal', @isnumeric);
addRequired(p, 'input_freq_hz', @(x) isnumeric(x) && x > 0);
addParameter(p, 'emg_fs', 2000, @(x) isnumeric(x) && x > 0);
addParameter(p, 'direction_axis', 4, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'cycle_type', 'full', @(x) ischar(x) || isstring(x));
parse(p, EMG_rms, imu_signal, input_freq_hz, varargin{:});

emg_fs     = p.Results.emg_fs;
rot_col    = p.Results.direction_axis;
cycle_type = lower(char(p.Results.cycle_type));

%% --- Orient Input Matrices ---
if size(EMG_rms, 1) ~= 64 && size(EMG_rms, 2) == 64
    EMG_rms = EMG_rms';
end
if size(imu_signal, 1) < size(imu_signal, 2)
    imu_signal = imu_signal';
end
if size(imu_signal, 2) >= rot_col
    imu_vec = imu_signal(:, rot_col);
else
    imu_vec = imu_signal(:, 1);
end

%% --- 1. Detect Debounced Zero-Crossings ---
imu_centered = imu_vec - mean(imu_vec);
raw_zc = find(imu_centered(1:end-1) .* imu_centered(2:end) <= 0);

if isempty(raw_zc)
    error('No zero crossings detected in the provided IMU signal.');
end

median_half_period = median(diff(raw_zc));
debounce_samples   = max(1, round(0.25 * median_half_period));

crossings_imu = raw_zc(1);
for k = 2:length(raw_zc)
    if (raw_zc(k) - crossings_imu(end)) >= debounce_samples
        crossings_imu = [crossings_imu; raw_zc(k)]; %#ok<AGROW>
    end
end

%% --- 2. Map Crossings to EMG & Form Cycle Boundaries ---
N_imu = length(imu_vec);
N_emg = size(EMG_rms, 2);
crossings_emg = 1 + round((crossings_imu - 1) * ((N_emg - 1) / (N_imu - 1)));

if strcmp(cycle_type, 'full')
    valid_starts = 1:2:(length(crossings_emg) - 2);
    cycle_starts = crossings_emg(valid_starts);
    cycle_ends   = crossings_emg(valid_starts + 2);
else % 'half'
    cycle_starts = crossings_emg(1 : end-1);
    cycle_ends   = crossings_emg(2 : end);
end

num_cycles = length(cycle_starts);
if num_cycles < 1
    error('Not enough zero-crossings to form at least one complete cycle.');
end

%% --- 3. Identify Outlier Cycles (>5% Deviation from Expected Period) ---
raw_cycle_lengths = (cycle_ends - cycle_starts) + 1; % Total samples per cycle

% Expected period parameters
expected_period_sec     = 1.0 / input_freq_hz;
expected_period_samples = expected_period_sec * emg_fs;
if strcmp(cycle_type, 'half')
    expected_period_samples = expected_period_samples / 2;
end

% 5% threshold allowance (in samples)
max_allowed_dev_samples = 0.05 * expected_period_samples;

% Evaluate deviation relative to expected period
cycle_deviations = abs(raw_cycle_lengths - expected_period_samples);
outlier_cycles   = find(cycle_deviations > max_allowed_dev_samples);
valid_cycles     = find(cycle_deviations <= max_allowed_dev_samples);

% Print Cycle Report to Command Window
fprintf('\n================ CYCLE ANALYSIS REPORT ================\n');
fprintf('Total Cycles Detected           : %d\n', num_cycles);
fprintf('Expected Period                 : %.4f s (%d samples at %d Hz)\n', ...
    expected_period_sec, round(expected_period_samples), emg_fs);
fprintf('5%% Allowed Deviation Threshold  : ±%.1f samples (%.4f s)\n', ...
    max_allowed_dev_samples, max_allowed_dev_samples / emg_fs);
fprintf('Valid Cycles Retained           : %d\n', length(valid_cycles));
fprintf('Outlier Cycles Flagged/Excluded : %d', length(outlier_cycles));
if ~isempty(outlier_cycles)
    fprintf(' (Cycle Indices: %s)\n', mat2str(outlier_cycles));
else
    fprintf('\n');
end
fprintf('=======================================================\n\n');

if isempty(valid_cycles)
    % Estimate actual frequency from detected cycle lengths (before erroring out)
    est_period_samples = median(raw_cycle_lengths);
    est_period_sec      = est_period_samples / emg_fs;
    est_freq_hz          = 1 / est_period_sec;
    if strcmp(cycle_type, 'half')
        est_freq_hz = est_freq_hz / 2; % half-cycles are twice the frequency of full cycles
    end

    fprintf('\n*** ALL CYCLES FLAGGED AS OUTLIERS ***\n');
    fprintf('Input frequency was %.4f Hz, but detected cycles suggest a period of\n', input_freq_hz);
    fprintf('%.1f samples (%.4f s), corresponding to an estimated frequency of %.4f Hz.\n', ...
        est_period_samples, est_period_sec, est_freq_hz);
    fprintf('Consider re-running with ''input_freq_hz'' ≈ %.4f.\n\n', est_freq_hz);
    disp(expected_period_samples)
    fprintf('actual')
    disp(raw_cycle_lengths)

    error('All cycles were flagged as outliers. Check your input frequency or signal quality. Estimated frequency: %.4f Hz.', est_freq_hz);
end

%% --- 4. Target Median Length Across Valid Cycles ---
median_len = round(median(raw_cycle_lengths(valid_cycles)));

%% --- 5. Segment & Interpolate Valid Cycles to L_median ---
n_valid = length(valid_cycles);
rms_cycles_cell = cell(64, n_valid);
target_grid     = linspace(0, 1, median_len);

for idx = 1:n_valid
    c = valid_cycles(idx);
    i_start = cycle_starts(c);
    i_end   = cycle_ends(c);
    
    raw_slice     = EMG_rms(:, i_start:i_end);
    n_raw_samples = size(raw_slice, 2);
    orig_grid     = linspace(0, 1, n_raw_samples);
    
    for ch = 1:64
        rms_cycles_cell{ch, idx} = interp1(orig_grid, raw_slice(ch, :), target_grid, 'pchip');
    end
end

%% --- 6. Compute Channel Means Across Valid Cycles ---
channel_mean_rms = zeros(64, median_len);
for ch = 1:64
    % Curly braces {ch, :} correctly generate a comma-separated list of vectors
    ch_cycles_matrix = vertcat(rms_cycles_cell{ch, :}); % [N_valid x L_median]
    channel_mean_rms(ch, :) = mean(ch_cycles_matrix, 1);
end

end