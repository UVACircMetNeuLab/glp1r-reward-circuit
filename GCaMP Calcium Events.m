clear all;
clc;

%% Parameters
samplerate = 120; % Hz
exclude_time = 5; % seconds to exclude from start
event_duration_threshold = 1.5; % seconds
event_samples_threshold = event_duration_threshold * samplerate;

%% === Load VEHICLE File ===
[veh_filename, veh_path] = uigetfile('*.doric', 'Select the VEHICLE recording file');
if isequal(veh_filename, 0)
    disp('User cancelled');
    return;
end
veh_data = ExtractDataAcquisition(fullfile(veh_path, veh_filename));

veh_time = veh_data(3).Data(2).Data;
veh_405 = veh_data(3).Data(1).Data;
veh_465 = veh_data(4).Data(1).Data;

veh_time = veh_time((exclude_time * samplerate + 1):end);
veh_405 = veh_405((exclude_time * samplerate + 1):end);
veh_465 = veh_465((exclude_time * samplerate + 1):end);

veh_fit = polyfit(veh_405, veh_465, 1);
veh_fitted_405 = polyval(veh_fit, veh_405);
veh_dff = ((veh_465 - veh_fitted_405) ./ veh_fitted_405) * 100;

veh_mean = mean(veh_dff);
veh_std = std(veh_dff);
zscore_vehicle = (veh_dff - veh_mean) / veh_std;
total_auc_veh = trapz(veh_time, zscore_vehicle);

% Event detection for vehicle
thresh_veh = max(median(zscore_vehicle) + 2 * std(zscore_vehicle), 1.5);
mask_veh = zscore_vehicle > thresh_veh;
starts_veh = strfind([0 mask_veh'], [0 1]);
ends_veh = strfind([mask_veh' 0], [1 0]);

durations_veh = ends_veh - starts_veh;
valid_veh = durations_veh >= event_samples_threshold;
starts_veh = starts_veh(valid_veh);
ends_veh = ends_veh(valid_veh);
event_times_veh = veh_time(starts_veh);

fprintf('\n--- VEHICLE Analysis (z-scored to self) ---\n');
fprintf('Vehicle AUC: %.2f\n', total_auc_veh);
fprintf('Vehicle Events: %d\n', length(starts_veh));
fprintf('Start times: %.2f, ', event_times_veh); fprintf('\n');

%% === Load DRUG File or Reuse VEHICLE ===
reuse = input('Use the same file for DRUG analysis? (1 = yes, 0 = load different): ');
if reuse == 1
    drug_time = veh_time;
    drug_405 = veh_405;
    drug_465 = veh_465;
else
    [drug_filename, drug_path] = uigetfile('*.doric', 'Select the DRUG recording file');
    if isequal(drug_filename, 0)
        disp('User cancelled');
        return;
    end
    drug_data = ExtractDataAcquisition(fullfile(drug_path, drug_filename));

    drug_time = drug_data(3).Data(2).Data;
    drug_405 = drug_data(3).Data(1).Data;
    drug_465 = drug_data(4).Data(1).Data;

    drug_time = drug_time((exclude_time * samplerate + 1):end);
    drug_405 = drug_405((exclude_time * samplerate + 1):end);
    drug_465 = drug_465((exclude_time * samplerate + 1):end);
end

% ΔF/F and z-score drug
drug_fit = polyfit(drug_405, drug_465, 1);
drug_fitted_405 = polyval(drug_fit, drug_405);
drug_dff = ((drug_465 - drug_fitted_405) ./ drug_fitted_405) * 100;
zscore_drug = (drug_dff - veh_mean) / veh_std;
total_auc_drug = trapz(drug_time, zscore_drug);

% Event detection for drug
thresh = max(median(zscore_drug) + 2 * std(zscore_drug), 1.5);
mask = zscore_drug > thresh;
starts = strfind([0 mask'], [0 1]);
ends = strfind([mask' 0], [1 0]);

durations = ends - starts;
valid = durations >= event_samples_threshold;
starts = starts(valid);
ends = ends(valid);
event_times = drug_time(starts);

fprintf('\n--- DRUG Analysis (z-scored to VEHICLE) ---\n');
fprintf('Drug AUC: %.2f\n', total_auc_drug);
fprintf('Drug Events: %d\n', length(starts));
fprintf('Start times: %.2f, ', event_times); fprintf('\n');

%% === Plot DRUG Trace with Events ===
figure;
subplot(3,1,1);
plot(drug_time, drug_405, 'b', drug_time, drug_465, 'r');
xlabel('Time (s)'); ylabel('Signal'); title('Raw Signals - DRUG'); legend('405', '465');

subplot(3,1,2);
plot(drug_time, zscore_drug, 'k'); xlabel('Time (s)'); ylabel('Z-score'); title('Z-scored ΔF/F - DRUG');

subplot(3,1,3);
plot(drug_time, zscore_drug, 'k'); hold on;
for i = 1:length(starts)
    idx = starts(i):ends(i);
    plot(drug_time(idx), zscore_drug(idx), 'r', 'LineWidth', 2);
end
plot(event_times, zscore_drug(starts), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
xlabel('Time (s)'); ylabel('Z-score'); title('Detected Events (Z > 2, ≥1.5s)');
hold off;

%% === Save to CSV for Heatmap (with first/append logic) ===
save_csv = input('Save z-scored DRUG trace to heatmap CSV? (1 = yes, 0 = no): ');
if save_csv == 1
    is_first = input('Is this the first mouse for this group? (1 = yes, 0 = no): ');
    
    if is_first == 1
        % Ask for folder and file name to create new CSV
        save_folder = uigetdir(pwd, 'Select folder to save new heatmap CSV');
        if isequal(save_folder, 0)
            disp('User cancelled.');
            return;
        end
        new_name = input('Enter new CSV file name (without .csv): ', 's');
        full_csv_path = fullfile(save_folder, [new_name, '.csv']);
        
        % Write new file with time + current z-score
        data_to_save = [drug_time(:), zscore_drug(:)];
        writematrix(data_to_save, full_csv_path);
        fprintf('New heatmap file created: %s\n', full_csv_path);
        
    else
        % Append to existing file
        [csv_file, csv_path] = uigetfile('*.csv', 'Select existing CSV to append');
        if isequal(csv_file, 0)
            disp('User cancelled.');
            return;
        end
        full_csv_path = fullfile(csv_path, csv_file);
        existing_data = readmatrix(full_csv_path);

        if size(existing_data, 1) ~= length(drug_time)
            error('Time mismatch between trace and existing file.');
        end

        % Append new z-score column
        updated_data = [existing_data, zscore_drug(:)];
        writematrix(updated_data, full_csv_path);
        fprintf('Z-score column appended to %s\n', full_csv_path);
    end
end



%% Select CSV files
[file_list, path] = uigetfile('*.csv', 'Select z-score CSV files to combine', 'MultiSelect', 'on');
if isequal(file_list, 0)
    disp('User cancelled');
    return;
end

if ischar(file_list)
    file_list = {file_list}; % convert to cell if only one file selected
end

num_files = length(file_list);
z_columns = cell(num_files, 1);
lengths = zeros(num_files, 1);

% First, read all files and store their z-score column and lengths
for i = 1:num_files
    full_path = fullfile(path, file_list{i});
    data = readmatrix(full_path);
    
    if size(data,2) < 2
        error('File %s has fewer than 2 columns. Expected: time + z-score.', file_list{i});
    end
    
    z = data(:,2);
    z_columns{i} = z;
    lengths(i) = length(z);
end

% Find the minimum trace length
min_length = min(lengths);

% Truncate and combine
z_matrix = zeros(min_length, num_files);
for i = 1:num_files
    z_matrix(:,i) = z_columns{i}(1:min_length);
end

%% Save combined matrix
[save_file, save_path] = uiputfile('combined_zscores.csv', 'Save combined z-score matrix as');
if isequal(save_file, 0)
    disp('User cancelled saving.');
else
    writematrix(z_matrix, fullfile(save_path, save_file));
    fprintf('Combined z-score matrix saved to: %s\n', fullfile(save_path, save_file));
end

%%
clear all;
clc;

%% Select combined z-score CSV
[file, path] = uigetfile('*.csv', 'Select combined z-score matrix');
if isequal(file, 0)
    disp('User cancelled');
    return;
end

full_path = fullfile(path, file);
z_matrix = readmatrix(full_path); % Each column = one mouse

% Transpose so: each row = one mouse
z_matrix_t = z_matrix';

%% Plot Heatmap
figure;
imagesc(z_matrix_t); % rows = mice, cols = time
colormap jet;
colorbar;

xlabel('Time (samples)');
ylabel('Mouse #');
title('Z-scored Calcium Activity Heatmap');

% Optional: add ticks
yticks(1:size(z_matrix_t, 1));
yticklabels(1:size(z_matrix_t, 1)); % Mouse numbers

figure;
imagesc(z_matrix_t); % rows = mice, cols = time
colormap jet;
colorbar;

xlabel('Time (samples)');
ylabel('Mouse #');
title('Z-scored Calcium Activity Heatmap');

% Add color scale limits
caxis([-4 4]); % <-- Set min and max of color scale

% Optional mouse number labels
yticks(1:size(z_matrix_t, 1));
yticklabels(1:size(z_matrix_t, 1));

%%clear all; clc;

%% === Select Files ===
[veh_file, veh_path] = uigetfile('*.csv', 'Select VEHICLE z-score CSV');
if isequal(veh_file, 0)
    error('Vehicle file not selected.');
end

[drug_file, drug_path] = uigetfile('*.csv', 'Select DRUG z-score CSV');
if isequal(drug_file, 0)
    error('Drug file not selected.');
end

vehicle_data = readmatrix(fullfile(veh_path, veh_file)); % rows = time, cols = mice
drug_data = readmatrix(fullfile(drug_path, drug_file));  % same structure

%% — Sanity-check & auto-align sizes — 
% vehicle_data and drug_data just read in

[rtV, cV] = size(vehicle_data);
[rtD, cD] = size(drug_data);
fprintf('Vehicle is %d×%d,  Drug is %d×%d\n', rtV, cV, rtD, cD);

% 1) If you have a leading time column, drop it:
if any(vehicle_data(:,1) == (0:rtV-1)') && any(drug_data(:,1) == (0:rtD-1)')
    vehicle_data = vehicle_data(:,2:end);
    drug_data    = drug_data(:,2:end);
    [rtV, cV] = size(vehicle_data);
    [rtD, cD] = size(drug_data);
    fprintf('Dropped time column → Vehicle %d×%d, Drug %d×%d\n', rtV, cV, rtD, cD);
end

% 2) Truncate to the smallest common size  
minRows = min(rtV, rtD);
minCols = min(cV, cD);

if rtV ~= rtD || cV ~= cD
    warning('Mismatched sizes — truncating to %d×%d', minRows, minCols);
    vehicle_data = vehicle_data(1:minRows, 1:minCols);
    drug_data    = drug_data(1:minRows, 1:minCols);
end

% now proceed with
[num_timepoints, num_mice] = size(vehicle_data);


[num_timepoints, num_mice] = size(vehicle_data);

if size(drug_data, 1) ~= num_timepoints || size(drug_data, 2) ~= num_mice
    error('Vehicle and Drug files do not match in size.');
end

%% === Parameters ===
smoothing_window = 10; % for movmean
n_low = 360;  % number of lowest points for min
n_high = 360; % number of highest points for max

%% === Normalize and Smooth ===
vehicle_norm = zeros(num_mice, num_timepoints);
drug_norm = zeros(num_mice, num_timepoints);

for i = 1:num_mice
    veh_trace = vehicle_data(:,i);

    % Compute min and max based on sorted vehicle z-score trace
    sorted_veh = sort(veh_trace);
    min_val = mean(sorted_veh(1:n_low));
    max_val = mean(sorted_veh(end-n_high+1:end));

    % Prevent divide-by-zero
    if max_val == min_val
        warning('Mouse %d: max == min in vehicle trace. Skipping this mouse.', i);
        continue;
    end

    % Normalize both VEHICLE and DRUG to vehicle-derived min/max
    veh_norm = (veh_trace - min_val) / (max_val - min_val);
    drug_norm_single = (drug_data(:,i) - min_val) / (max_val - min_val);

    % Smooth with moving average
    veh_smoothed = movmean(veh_norm, smoothing_window);
    drug_smoothed = movmean(drug_norm_single, smoothing_window);

    % Store in matrix (each row = one mouse)
    vehicle_norm(i,:) = veh_smoothed;
    drug_norm(i,:) = drug_smoothed;
end


%% === Plot Heatmaps ===
figure;

% --- Vehicle panel ---
ax1 = subplot(2,1,1);
imagesc(vehicle_norm);        % draw the vehicle heatmap
set(ax1, 'YDir', 'normal');   % orient rows top→bottom
axis tight;                   % auto-tighten the axes
colormap(ax1, 'jet');         % use the jet colormap
caxis([0 1]);                 % keep your 0–1 scale
colorbar;                     
title('Vehicle (Min–Max Normalized per Mouse)');
xlabel('Timepoints');
ylabel('Mouse');

% --- Drug panel ---
ax2 = subplot(2,1,2);
imagesc(drug_norm);           
set(ax2, 'YDir', 'normal');
axis tight;
colormap(ax2, 'jet');
caxis([0 1]);
colorbar;
title('Drug (Normalized using Vehicle Range)');
xlabel('Timepoints');
ylabel('Mouse');



%% Select combined z-score CSV
[orig_file, orig_path] = uigetfile('*.csv', 'Select combined z-score CSV');
if isequal(orig_file,0)
    disp('User cancelled');
    return;
end
orig_full = fullfile(orig_path, orig_file);

% Read in your full-resolution matrix (rows = time, cols = mice)
z_matrix = readmatrix(orig_full);

%% Downsampling parameters
% Prompt for factor or hard-code it
dsFactor = input('Enter downsample factor (e.g. 10 to keep 1/10th of points): ');

% Simple pick-every-Nth approach:
z_ds = z_matrix(1:dsFactor:end, :);

% Alternatively, for anti-aliased downsampling:
% z_ds = zeros(floor(size(z_matrix,1)/dsFactor), size(z_matrix,2));
% for i = 1:size(z_matrix,2)
%     z_ds(:,i) = decimate(z_matrix(:,i), dsFactor);
% end

%% Save the downsampled matrix
[save_file, save_path] = uiputfile('downsampled_zscores.csv', ...
                                   'Save downsampled matrix as');
if isequal(save_file,0)
    disp('User cancelled saving.');
else
    writematrix(z_ds, fullfile(save_path, save_file));
    fprintf('Downsampled matrix saved to:\n%s\n', fullfile(save_path, save_file));
end