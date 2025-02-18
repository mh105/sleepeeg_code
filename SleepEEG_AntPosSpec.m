function [] = SleepEEG_AntPosSpec(subID, fnsuffix, project, EEG, chanlist, fileID, outputDir)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - ANTPOSSPEC**
%
% - used to plot the time traces of an anterior lead (84-Z2) and a
% posterior lead (25-Z10) and the multi-taper spectrograms. Outputs are
% saved to the same OutputDir
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
%           - chanlist:     Channel numbers to compute the spectrograms for.
%                           default: [25, 84]
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
disp('-------------------------------')
disp([mHead 'SleepEEG_AntPosSpec()']);
disp('-------------------------------')

%% Load Data if the EEG structure is not an input
if nargin < 4
    [ dataDir, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
    assert(isfile(fullfile(dataDir, subID, 'set', [subID, '_', fnsuffix, '_ds500_Z3.set'])),...
        [mHead, 'Downsampled data set .set file is not available for ' fileID])
    EEG = SleepEEG_loadset(subID, fnsuffix, project);
end

if ~exist('chanlist', 'var')
    % By default, this is going to be 25-Z10 (posterior lead) and 84-Z2
    % (anterior lead) for visualizing the resting state recordings
    chanlist = [25, 84];
end

if ~exist('fileID', 'var')
    [ ~, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

fileID = [fileID, '_', EEG.refscheme]; % specify the referencing scheme in filenames

% Report plotted channels
disp([mHead, 'Plotting channels: ', num2str(chanlist)])

%% Prepare for plotting and labelling
try
    %     EEG.endimp(SleepEEG_impcheck(subID, fnsuffix, EEG, fileID, outputDir, false)) = NaN;
    implist = EEG.endimp(chanlist); % grab the end impedance value
catch
    implist = NaN*zeros(size(chanlist));
end

% check whether unipolar reference is being used
ref_channel = [];
for ii = 1:length(EEG.chanlocs)
    if strcmp(EEG.refscheme, EEG.chanlocs(ii).labels)
        ref_channel = ii;
    end
end
if ~isempty(ref_channel)
    chanlist = [chanlist, ref_channel];
else
    chanlist = [chanlist, 87];  % just to place an electrode on topoplot
end

%% Visualize sensor positions
figure; topoplot([],EEG.chanlocs(chanlist),'style','both','electrodes','ptslabels','emarker', {'.', 'k', 15, 1});
L = findobj(gcf, 'type', 'Text');
for ind = 1:length(chanlist)-1
    set(L(length(chanlist)+1-ind), 'FontSize', 20)
    set(L(length(chanlist)+1-ind), 'Color', [0 0 0])
end
% indicate the reference scheme
if ~isempty(ref_channel)
    set(L(1), 'FontSize', 20)
    set(L(1), 'Color', [1,0,0])
    set(L(1), 'String', [L(1).String, '(Ref)'])
else
    set(L(1), 'FontSize', 20)
    set(L(1), 'Color', [1,0,0])
    set(L(1), 'String', [EEG.refscheme, ' (Ref)'])
end
% figure specifications
title(['Channel locations for ' EEG.chanlocs(chanlist(1)).labels '-' num2str(chanlist(1))...
    ', ' EEG.chanlocs(chanlist(2)).labels '-' num2str(chanlist(2))], 'Interpreter', 'none')
legend([EEG.chanlocs(chanlist(1)).labels ' IMP= ' num2str(implist(1))],...
    [EEG.chanlocs(chanlist(2)).labels ' IMP= ' num2str(implist(2))])
legend('Location', 'best')
set(gca,'FontSize', 16)
saveas(gcf, fullfile(outputDir, [fileID '_two_channels_location.png']))

chanlist(end) = []; % remove the reference channel after creating topoplot

%% Visualize detrended electrode time traces
t = EEG.times./(1000); % convert to second scale
figure
for i = chanlist
    plot(t, detrend(EEG.data(i,:)))
    hold on
end
hold off
title(['Detrended Trace for Sensor ' EEG.chanlocs(chanlist(1)).labels '-' num2str(chanlist(1))...
    ', ' EEG.chanlocs(chanlist(2)).labels '-' num2str(chanlist(2))], 'Interpreter', 'none')
legend([EEG.chanlocs(chanlist(1)).labels ' IMP= ' num2str(implist(1))],...
    [EEG.chanlocs(chanlist(2)).labels ' IMP= ' num2str(implist(2))])
legend('Location', 'best')
xlabel('Time(sec)')
ylabel('\mu V')
set(gca, 'FontSize', 16)
saveas(gcf, fullfile(outputDir, [fileID '_two_channels_detrend_trace.png']))

%% Compute spectrograms of the two channels

mt_spectrogram = cell(2,length(chanlist));
stimes = cell(2,length(chanlist));
sfreqs = cell(2,length(chanlist));

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

% compute spectrogram up to Nyquist freq
for i = 1:length(chanlist)
    % display MTM parameters only once
    if i == 1; verbose=true; else; verbose=false; end

    [mt_spectrogram{1,i},stimes{1,i},sfreqs{1,i}] = multitaper_spectrogram_mex(EEG.data(chanlist(i),:), EEG.srate, [0, (EEG.srate)/2],...
        [5, 9], [5, 1], 0, 'linear', 'unity', false, verbose);
end

% Repeat spectral estimation after applying a band-stop filter around 60Hz
if EEG.srate / 2 > 55
    for i = 1:length(chanlist)
        eeg = bandstop(EEG.data(chanlist(i),:), [55, min(65, EEG.srate / 2)], EEG.srate);
        [mt_spectrogram{2,i},stimes{2,i},sfreqs{2,i}] = multitaper_spectrogram_mex(eeg, EEG.srate, [0, (EEG.srate)/2],...
            [5, 9], [5, 1], 0, 'linear', 'unity', false, false);
    end
else
    mt_spectrogram(2,:) = mt_spectrogram(1,:);
    stimes(2,:) = stimes(1,:);
    sfreqs(2,:) = sfreqs(1,:);
end

%% Plot MTM spectrograms
% Generate a plot
figure
ax = figdesign(4,3,'merge',{1:2, 4:5, 7:8, 10:11});
for ii = 1:length(ax); title(ax(ii), ii); end
set(gcf, 'units', 'pixels', 'Position', [0 0 1400 1000]);

% For the spectrogram up to Nyquist, don't actually plot it all the way to
% Nyquist because low-pass filter was applied before downsampling such that
% there will be a noise drop off beyond 220Hz. This will mess up the color
% scaling on the spectrograms. So we will just plot it up to 210Hz.
fullspecfreq = [0, 210];
fullspecfreq_idx = sfreqs{1,1}>=fullspecfreq(1) & sfreqs{1,1}<=fullspecfreq(2);

% For zooming in on frequency 1:40Hz
visfreq = [0, 40];
freq_idx = sfreqs{2,1}>=visfreq(1) & sfreqs{2,1}<=visfreq(2);


% Plot the 1st electrode
axes(ax(1));
imagesc(stimes{1,1},sfreqs{1,1}(fullspecfreq_idx),pow2db(mt_spectrogram{1,1}(fullspecfreq_idx,:)));
axis xy
% xlabel('Time (s)');
ylabel('Frequency (Hz)');
c=colorbar;
title(['Spectrogram ' EEG.chanlocs(chanlist(1)).labels '-' num2str(chanlist(1)) ' posterior'])
ylabel(c,'PSD (dB)');
set(gca, 'FontSize', 16)
colormap jet
climscale;
axis tight
cspect1 = clim;


axes(ax(2));
y = median(mt_spectrogram{1,1}, 2);
plot(sfreqs{1,1}, pow2db(y), 'LineWidth', 2)
hold on
legend([EEG.chanlocs(chanlist(1)).labels ' IMP= ' num2str(implist(1))])
legend('Location', 'best')
title(['Median Spectrum ' EEG.chanlocs(chanlist(1)).labels '-' num2str(chanlist(1))])
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB)')
set(gca, 'FontSize', 16)
axis tight
xlim(fullspecfreq)
ylimit1 = ylim;


% Spectrogram zoomed in on frequency 1:40Hz
axes(ax(3));
imagesc(stimes{2,1},sfreqs{2,1}(freq_idx),pow2db(mt_spectrogram{2,1}(freq_idx,:)));
axis xy
% xlabel('Time (s)');
ylabel('Frequency (Hz)');
c=colorbar;
title(['Spectrogram ' EEG.chanlocs(chanlist(1)).labels '-' num2str(chanlist(1))])
ylabel(c,'PSD (dB)');
set(gca, 'FontSize', 16)
colormap jet
climscale;
axis tight


% averaged spectra 0-40Hz for both electrodes overlaid on the same plot
axes(ax(4));
y = median(mt_spectrogram{2,1}, 2);
plot(sfreqs{2,1}, pow2db(y), 'LineWidth', 2)
hold on
y = median(mt_spectrogram{2,2}, 2);
plot(sfreqs{2,2}, pow2db(y), 'LineWidth', 2)
legend([EEG.chanlocs(chanlist(1)).labels ' IMP= ' num2str(implist(1))],...
    [EEG.chanlocs(chanlist(2)).labels ' IMP= ' num2str(implist(2))])
legend('Location', 'best')
title('Median Spectra (0-40 Hz)')
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB)')
set(gca, 'FontSize', 16)
axis tight
xlim(visfreq)


% Plot the 2nd electrode
axes(ax(5));
imagesc(stimes{1,2},sfreqs{1,2}(fullspecfreq_idx),pow2db(mt_spectrogram{1,2}(fullspecfreq_idx,:)));
axis xy
% xlabel('Time (s)');
ylabel('Frequency (Hz)');
c=colorbar;
title(['Spectrogram ' EEG.chanlocs(chanlist(2)).labels '-' num2str(chanlist(2)) ' anterior'])
ylabel(c,'PSD (dB)');
set(gca, 'FontSize', 16)
colormap jet
climscale;
axis tight
cspect2 = clim;


axes(ax(6));
y = median(mt_spectrogram{1,2}, 2);
plot(sfreqs{1,2}, pow2db(y), 'LineWidth', 2)
hold on
legend([EEG.chanlocs(chanlist(2)).labels ' IMP= ' num2str(implist(2))])
legend('Location', 'best')
title(['Median Spectrum ' EEG.chanlocs(chanlist(2)).labels '-' num2str(chanlist(2))])
xlabel('Frequency (Hz)')
% ylabel('PSD (dB])')
set(gca, 'FontSize', 16)
axis tight
xlim(fullspecfreq)
ylimit2 = ylim;


% Spectrogram zoomed in on frequency 1:40Hz
axes(ax(7));
imagesc(stimes{2,2},sfreqs{2,2}(freq_idx),pow2db(mt_spectrogram{2,2}(freq_idx,:)));
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
c=colorbar;
title(['Spectrogram ' EEG.chanlocs(chanlist(2)).labels '-' num2str(chanlist(2))])
ylabel(c,'PSD (dB)');
set(gca, 'FontSize', 16)
colormap jet
climscale;
axis tight


% Information panel
axes(ax(8));
title('Notes', 'FontSize', 16)
line1 = ['Task: ' subID, '_', fnsuffix];
line1 = strrep(line1,'_','-');
text(0.1, 0.8, line1, 'FontSize', 16)
line2 = ['Channels: ' EEG.chanlocs(chanlist(1)).labels,'-',EEG.refscheme, ', ',...
    EEG.chanlocs(chanlist(2)).labels,'-',EEG.refscheme];
text(0.1, 0.5, line2, 'FontSize', 16)
line3 = 'MTM params: [5,9], [5,1]';
text(0.1, 0.2, line3, 'FontSize', 16)


% Adjust the plot axes limits
axes(ax(1));
clim([min(cspect1(1), cspect2(1)), max(cspect1(2), cspect2(2))]);
axes(ax(5));
clim([min(cspect1(1), cspect2(1)), max(cspect1(2), cspect2(2))]);
axes(ax(2));
ylim([min(ylimit1(1), ylimit2(1)), max(ylimit1(2), ylimit2(2))]);
axes(ax(6));
ylim([min(ylimit1(1), ylimit2(1)), max(ylimit1(2), ylimit2(2))]);


% Link axes
linkaxes([ax(1), ax(3), ax(5), ax(7)], 'x')

saveas(gcf, fullfile(outputDir, [fileID '_two_channels_spectrogram.png']))
% Also save a MATLAB .fig file for manipulations and zoom-in/out
savefig(fullfile(outputDir, [fileID '_two_channels_spectrogram.fig']))

close all

end
