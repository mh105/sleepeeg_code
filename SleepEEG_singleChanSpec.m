function [] = SleepEEG_singleChanSpec(subID, fnsuffix, data, channelNum, project, srate, hrscale, fileID, outputDir)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - SINGLECHANSPEC**
%
% - used to plot the time traces and spectrograms of several channels
% without downsampling. Typically such data are the outputs of functions
% SleepEEG_extractChan.m followed by SleepEEG_reconChan.m
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep".
%
%                           ***the final file name is in the form: subID_fnsuffix.cnt***
%
%           - data:         a M x N matrix containing data of one or more
%                           channels without downsampling (could be
%                           downsampled as well).
%
%           - channelNum:   a 1 x M vector containing channel numbers.
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
%           - srate:        original sampling rate, double type.
%                           default: 1000 (Hz)
%
%           - hrscale:      whether to use hour scale on time axis, default
%                           is false, which uses seconds.
%                           default: false
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

% Check the dimensions of data and channelNum agree
assert(size(data,1) == length(channelNum), 'Number of channels in data dismatches the number in channelNum!')

% Default flag variables
if nargin < 5
    project = '';
    srate = 1000; % default sampling rate is 1000Hz
    hrscale = false;
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';

%%
disp('-----------------------------------')
disp([mHead 'SleepEEG_singleChanSpec()']);
disp('-----------------------------------')

%% Report extracted channels
disp([mHead, 'Single channels plotted in original sampling frequency: ', num2str(channelNum)])

%% Define Directories of Codes and Data Folders
if ~exist('fileID', 'var')
    [ ~, ~, fileID ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

if ~exist('outputDir', 'var')
    [ ~, ~, ~, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

%% Plot single channel spectrograms
for i = 1:length(channelNum) % loop through each channel number

    %% Compute spectrograms of the channel
    %multitaper_spectrogram_mex.m
    %   Input:
    %   data: <number of samples> x 1 vector - time series data -- required
    %   Fs: double - sampling frequency in Hz  -- required
    %   frequency_range: 1x2 vector - [<min frequency>, <max frequency>] (default: [0 nyquist])
    %   taper_params: 1x2 vector - [<time-halfbandwidth product>, <number of tapers>] (default: [5 9])
    %   window_params: 1x2 vector - [window size (seconds), step size (seconds)] (default: [30 5])
    %   min_NFFT: double - minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x) (default: 0)
    %   detrend_opt: string - detrend data window ('linear' (default), 'constant', 'off');
    %   weighting: string - weighting of tapers ('unity' (default), 'eigen', 'adapt');
    %   plot_on: boolean to plot results (default: true)
    %   verbose: boolean to display spectrogram properties (default: true)

    % display MTM parameters only once
    if i == 1; verbose=true; else; verbose=false; end

    % compute spectrogram up to Nyquist freq
    [mt_spectrogram,stimes,sfreqs] = multitaper_spectrogram_mex(data(i,:), srate, [0, srate/2],...
        [5, 9], [30, 5], 0, 'constant', 'unity', false, verbose);

    %% Plot MTM spectrograms
    % Generate a plot
    figure
    ax = figdesign(3,2,'merge',{1:2, 3:4, 5:6});
    for ii = 1:length(ax); title(ax(ii), ii); end
    set(gcf, 'units', 'pixels', 'Position', [0 0 1400 1000]);


    axes(ax(1)); %#ok<*LAXES>
    imagesc(stimes, sfreqs, pow2db(mt_spectrogram));
    axis xy
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    c=colorbar;
    title(['Spectrogram for sensor ' num2str(channelNum(i))])
    ylabel(c,'PSD (dB)');
    set(gca, 'FontSize', 16)
    colormap jet
    climscale;
    axis tight

    %% Plot time trace on the same plot
    if hrscale
        t = linspace(0, length(data(i,:))/srate, length(data(i,:)))./60/60; % convert to hour scale
    else
        t = linspace(0, length(data(i,:))/srate, length(data(i,:))); % in second scale
    end

    % Visualize detrended electrodes
    axes(ax(2));
    plot(t, detrend(data(i,:)))
    hold on
    title(['Detrended Trace for sensor ' num2str(channelNum(i))])
    xlabel('Time (s)')
    ylabel('\mu V')
    set(gca, 'FontSize', 16)
    axis tight
    colorbar

    % Link axis to spectrogram
    linkaxes([ax(1) ax(2)], 'x')

    %% Plot the median spectrum
    axes(ax(3));
    y = median(mt_spectrogram, 2);
    plot(sfreqs, pow2db(y), 'LineWidth', 2)
    hold on
    title(['Median Spectrum for sensor ' num2str(channelNum(i))])
    xlabel('Frequency (Hz)')
    ylabel('Average PSD (dB)')
    set(gca, 'FontSize', 16)
    axis tight

    %% Save figure
    saveas(gcf, fullfile(outputDir, [fileID '_Raw_channel_' num2str(channelNum(i)) '_srate' num2str(srate) '_spectrogram.png']))
    % Also save a MATLAB .fig file for manipulations and zoom-in/out
    savefig(fullfile(outputDir, [fileID '_Raw_channel_' num2str(channelNum(i)) '_srate' num2str(srate) '_spectrogram.fig']))

    close all
end

end
