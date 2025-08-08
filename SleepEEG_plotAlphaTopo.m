function [] = SleepEEG_plotAlphaTopo(subID, fnsuffix, project, EEG, alpha_range, fileID, outputDir)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - PLOTALPHATOPO**
%
% - plot the topography of alpha activity during the recording
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep".
%
%                           ***the final file name is in the form: subID_fnsuffix.cnt***
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
%           - EEG:          data structure containing the EEG data and
%                           impedance values.
%
%           - alpha_range:  frequency range used to sum PSD for alpha activity.
%                           default: [8, 12]
%
%           - fileID:       name of the .cnt file (no .cnt suffix).
%
%           - outputDir:    directory path to folder for dumping all
%                           analyses outputs.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';

%%
disp('----------------------------------')
disp([mHead 'SleepEEG_plotAlphaTopo()']);
disp('----------------------------------')

%% Load Data if the EEG structure is not an input
if nargin < 4
    [ dataDir, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
    assert(isfile(fullfile(dataDir, subID, 'set', [subID, '_', fnsuffix, '_ds500_Z3.set'])),...
        [mHead, 'Downsampled data set .set file is not available for ' fileID])
    EEG = SleepEEG_loadset(subID, fnsuffix, project);
end

if ~exist('alpha_range', 'var') || isempty(alpha_range)
    alpha_range = [8, 12];
end

if ~exist('fileID', 'var')
    [ ~, ~, fileID ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

if ~exist('outputDir', 'var')
    [ ~, ~, ~, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

%% Check if the current recording is resting state
events = EEG.event;
task_start_code_idx = find(strcmp({events.type}, '100'));
task_code = events(task_start_code_idx(end) + 1).type;

if strcmp(task_code, '101')
    % take out the resting state eyes closed data
    % find the '2' trigger index
    trigger_2_idx = find(cellfun(@(x) strcmp(x, '2'), {EEG.event.type}), 1, 'last');
    index_2 = EEG.event(trigger_2_idx).latency;
    % '3' trigger index
    trigger_3_idx = find(cellfun(@(x) strcmp(x, '3'), {EEG.event.type}), 1, 'last');
    index_3 = EEG.event(trigger_3_idx).latency;
else
    index_2 = 1;
    index_3 = size(EEG.data, 2);
end

data = EEG.data(:, index_2:index_3);

%% MTM spectral analysis
alpha_power = zeros(1, size(data, 1));
for ii = 1:length(alpha_power)
    [mt_spectrogram,~,sfreqs] = multitaper_spectrogram_mex(data(ii, :), EEG.srate, [0, (EEG.srate)/2],...
        [5, 9], [5, 1], 0, 'linear', 'unity', false, false);
    mean_spectrum = mean(mt_spectrogram, 2);
    valid_sfreqs = sfreqs >= alpha_range(1) & sfreqs <= alpha_range(2);
    alpha_power(ii) = sum(mean_spectrum(valid_sfreqs));
end

%% Topoplot of alpha distribution across the scalp
f = figure;
subplot(1,2,1)
topoplot(alpha_power,EEG.chanlocs,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});
cb = colorbar;
cb.Label.String = '\muV^2/Hz';
cb.Label.FontSize = 16;
title({fileID, 'Alpha PSD Topoplot (Linear Scale)'}, 'Interpreter', 'none')
set(gca, 'FontSize', 16)
subplot(1,2,2)
topoplot(pow2db(alpha_power),EEG.chanlocs,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});
cb = colorbar;
cb.Label.String = 'dB';
cb.Label.FontSize = 16;
title({fileID, 'Alpha PSD Topoplot (Log Scale)'}, 'Interpreter', 'none')
set(gca, 'FontSize', 16)

saveas(gcf, fullfile(outputDir, [fileID '_Alpha_topoplot.png']))
close(f)

end
