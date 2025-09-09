clear all; clc;

%%This code has been used for GLP1R-dopamine experiments in the paper.
%% Define Parameters
samplerate = 120; % Sampling rate in Hz

% Prompt user to input time stamps for z-score calculation 
calc_time_stamps = input('Enter the list of time stamps for z-score calculation in seconds (e.g., [450, 900, 1350]): ');
plot_time_stamps = input('Enter the list of time stamps for z-score plotting in seconds (e.g., [500, 1000, 1500]): ');

% Exclude the first 5 seconds from all trial timestamps
calc_time_stamps = calc_time_stamps - 5;
plot_time_stamps = plot_time_stamps - 5;

% Check if any trial timestamps are negative after adjustment
if any(calc_time_stamps < 0) || any(plot_time_stamps < 0)
    error('Adjusted trial timestamps cannot be negative. Please check the video recording analysis.');
end

% Check if the inputs are non-empty numeric arrays
if isempty(calc_time_stamps) || ~isnumeric(calc_time_stamps) || isempty(plot_time_stamps) || ~isnumeric(plot_time_stamps)
    error('Invalid time stamps input. Please enter numeric arrays.');
end

% Prompt user for the time window around the plot time stamp
time_before_plot = input('Enter the time window before the plot time stamp in seconds: ');
time_after_plot = input('Enter the time window after the plot time stamp in seconds: ');

%% Load Data
[filename, path] = uigetfile('*.doric', 'Select the data file');
if isequal(filename, 0)
    disp('User selected Cancel');
    return;
end

% Load .doric file with error handling
try
    data = ExtractDataAcquisition(fullfile(path, filename));
catch
    error('Error loading the file. Please check the file format.');
end

% Extract signal data assuming 405 nm and 465 nm channels
try
    time_405nm = data(3).Data(2).Data;
    signal_405nm = data(3).Data(1).Data;
    time_465nm = data(4).Data(2).Data;
    signal_465nm = data(4).Data(1).Data;
catch
    error('Error extracting data channels. Please check the data structure.');
end

% Exclude the first 5 seconds of the recording
exclude_time = 5; % Time to exclude from the start in seconds
exclude_samples = round(exclude_time * samplerate);
time_405nm = time_405nm(exclude_samples:end);
signal_405nm = signal_405nm(exclude_samples:end);
time_465nm = time_465nm(exclude_samples:end);
signal_465nm = signal_465nm(exclude_samples:end);

% Use the original time and data
processed_time = time_405nm;
processed_signal_405nm = signal_405nm;
processed_signal_465nm = signal_465nm;

%% Fit Signals
% Fit 405 nm signals to 465 nm signals using linear least squares fit
fitted_405nm = polyfit(processed_signal_405nm, processed_signal_465nm, 1);
fitted_signal_405nm = polyval(fitted_405nm, processed_signal_405nm);

%% Calculate ΔF/F
delta_f_over_f = (processed_signal_465nm - fitted_signal_405nm) ./ fitted_signal_405nm;

%% Plot %ΔF/F
figure;
plot(processed_time, delta_f_over_f * 100, 'b');
xlabel('Time (s)');
ylabel('%ΔF/F');
title('%ΔF/F over Time');
grid on;



%% Plot Z-Scores around Plot Time Stamps
for i = 1:length(plot_time_stamps)
    % … your existing z-score plotting code …
end


%% Initialize Z-Score Array for Each Timestamp
z_scores_all = zeros(length(calc_time_stamps), length(delta_f_over_f)); % Preallocate array to store z-scores for each timestamp

% Initialize arrays to store baseline means and standard deviations
baseline_means = zeros(1, length(calc_time_stamps));
baseline_stds = zeros(1, length(calc_time_stamps));

%% Calculate Z-Scores for Each Calculation Time Stamp
for i = 1:length(calc_time_stamps)
    % Define the current calculation time stamp and the baseline window
    current_time_stamp = calc_time_stamps(i);
    baseline_start_time = current_time_stamp - 30; % 60 sec before trial start
    baseline_end_time = current_time_stamp - 0; % 30 sec window

    % Find the indices for the baseline period
    baseline_start_index = find(processed_time >= baseline_start_time, 1);
    baseline_end_index = find(processed_time >= baseline_end_time, 1);

    % Ensure the indices are valid
    if isempty(baseline_start_index) || isempty(baseline_end_index) || baseline_end_index <= baseline_start_index
        error('Invalid baseline time window for time stamp %d.', current_time_stamp);
    end

    % Calculate baseline mean and std
    baseline_means(i) = mean(delta_f_over_f(baseline_start_index:baseline_end_index));
    baseline_stds(i) = std(delta_f_over_f(baseline_start_index:baseline_end_index));

    % Calculate z-scores using the entire recording
    z_scores_all(i, :) = (delta_f_over_f - baseline_means(i)) / baseline_stds(i);
end


%% Prompt for Save Directory
save_path = uigetdir(path, 'Select the folder to save plots and data');
if save_path == 0
    disp('User selected Cancel');
    return;
end

%% Plot Z-Scores around Plot Time Stamps
for i = 1:length(plot_time_stamps)
    % Define the current plot time stamp
    current_plot_time_stamp = plot_time_stamps(i);

    % Find the indices for the plot time stamp
    plot_time_index = find(processed_time >= current_plot_time_stamp, 1);

    % Ensure the index is valid
    if isempty(plot_time_index)
        error('Invalid plot time stamp %d.', current_plot_time_stamp);
    end

    % Define the window around the plot time stamp
    window_start_index = max(1, plot_time_index - round(time_before_plot * samplerate));
    window_end_index = min(length(processed_time), plot_time_index + round(time_after_plot * samplerate));

    % Extract the z-scores and time for the window around the plot time stamp
    z_scores_window = z_scores_all(i, window_start_index:window_end_index);
    time_window = processed_time(window_start_index:window_end_index);

    % Plot z-scores with the plot time stamp at zero
    figure;
    plot(time_window - current_plot_time_stamp, z_scores_window, 'r');
    xlabel('Time relative to Plot Time Stamp (s)');
    ylabel('Z-Score');
    title(sprintf('Z-Scores around Plot Time Stamp %d', current_plot_time_stamp));
    grid on;

    % Save the plot
    saveas(gcf, fullfile(save_path, sprintf('z_scores_plot_time_stamp_%d.png', current_plot_time_stamp)));

    % Save time and z-score data into csv
    if i == 1
        % Initialize CSV data with time column
        csv_data = [time_window - current_plot_time_stamp, z_scores_window'];
    else
        % Append z-score data as new columns
        csv_data = [csv_data, z_scores_window'];
    end
end

% Calculate the average z-score across all trials
average_z_scores = mean(csv_data(:, 2:end), 2);

% Get the base filename (without extension) from the original data file
[~, base_filename, ~] = fileparts(filename);

% Append the average z-scores as the last column
csv_data = [csv_data, average_z_scores];

% Write CSV data to file using the same base name as the data file
csv_filename = fullfile(save_path, [base_filename, '_z_scores_trial_data.csv']);
csv_header = ['Time', arrayfun(@(x) sprintf('Trial_%d', x), 1:length(plot_time_stamps), 'UniformOutput', false), 'Average'];
csv_data_table = array2table(csv_data, 'VariableNames', csv_header);
writetable(csv_data_table, csv_filename);

disp(['Z-scores for each trial and average z-scores saved successfully as ', csv_filename]);

%% 

% Save Processed Data
processed_data_filename = fullfile(path, [filename(1:end-6), '_processed_data.mat']);
save(processed_data_filename, 'processed_time', 'delta_f_over_f', 'z_scores_all');
disp(['Processed data saved as ', processed_data_filename]);

% Display message
disp('Processed data, plots, and CSV saved successfully.');