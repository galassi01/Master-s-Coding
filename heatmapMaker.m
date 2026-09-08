% 1. Define dimensions
rows = 16;
cols = 6;
total_elements = rows * cols;

% 2. Create a medium-strength top-left bias using controlled decay
[X, Y] = meshgrid(1:cols, 1:rows);
distance_from_tl = sqrt((X - 1).^2 + (Y - 1).^2);
% Invert it so smaller distances have higher weights
tl_bias = 1 ./ (distance_from_tl + 1);

% 3. Generate random noise blended with our medium bias
rng('shuffle');
raw_noise = rand(rows, cols) .* (1 + 1.2 * medium_bias);

% Apply the spatial smoothing filter
smoothing_filter = [0.1, 0.2, 0.1; 0.2, 0.5, 0.2; 0.1, 0.2, 0.1];
smoothed_noise = conv2(raw_noise, smoothing_filter, 'same');

% 4. Sort the values to find exact tertile split points (1/3 each)
[~, sort_idx] = sort(smoothed_noise(:));

% Create the final matrix scaled between 0 and 1
heatmap_data = zeros(rows, cols);

% Group 1: Bottom 1/3 -> Map from 0 to 0.30
idx_low = sort_idx(1:32);
heatmap_data(idx_low) = linspace(0.20, 0.30, 32);

% Group 2: Middle 1/3 -> Map from 0.33 to 0.40
idx_mid = sort_idx(33:64);
heatmap_data(idx_mid) = linspace(0.33, 0.50, 32);

% Group 3: Top 1/3 (Threshold) -> Map from 0.45 to 1.00
idx_high = sort_idx(65:96);
heatmap_data(idx_high) = linspace(0.55, 1.00, 32);

% Shuffle locally within groups to keep a natural look
heatmap_data(idx_low) = heatmap_data(idx_low(randperm(32)));
heatmap_data(idx_mid) = heatmap_data(idx_mid(randperm(32)));
heatmap_data(idx_high) = heatmap_data(idx_high(randperm(32)));


% 5. Plotting the Heatmap
figure('Position', [200, 100, 450, 650]);
imagesc(heatmap_data);
colormap('parula'); 
colorbar;
caxis([0, 1]); 

% Formatting Axes
xticks(1:cols);
yticks(1:rows);
xlabel('Electrode');
ylabel('');
title('d Prime Motor Thresholds');
set(gca, 'YDir', 'reverse'); % Keep Row 1 at the top

% Note: Quadrant lines (xline/yline) have been completely removed