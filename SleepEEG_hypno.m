function [] = SleepEEG_hypno(subID, fnsuffix, channel, project, datasource)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - HYPNO**
%
% - this is a function used to extract the scored sleep stages from
% exported .txt files from the clinical Natus system and plot hypnograms of
% various sorts
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep".
%
%                           *** this function will use subID and fnsuffix
%                               to load the .set file in the set folder if
%                               'eeg' or 'both' is used as datasource:
%                               subID_fnsuffix_ds500_Z3.set ***
%
%           - channel:      a string character or a cell with different
%                           channels when datasource is 'both'. If
%                           datasource is 'eeg', then channel can be either
%                           a string or a double type channel number. If
%                           datasource is 'edf', then channel must be a
%                           string character from F3,F4,C3,C4,O1,O2.
%                           default: {'Z6', 'C3'}, or {87, 'C3'}
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
%           - datasource:   datasource to plot the spectrogram from. Can be
%                           'eeg' for grabbing channel from HD-EEG .set
%                           file, 'edf' for grabbing channel from aligned
%                           EDF+ file, or 'both' for plotting two
%                           spectrograms, one for each datasource.
%                           default: 'both'
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 4
    project = '';
    datasource = 'both';
elseif nargin < 5
    datasource = 'both';
end
if ~exist('channel', 'var')
    eegchannel = 'Z6';
    edfchannel = 'C3';
elseif isempty(channel)
    eegchannel = 'Z6';
    edfchannel = 'C3';
else
    if isscalar(channel) && strcmp(datasource, 'eeg')
        eegchannel = channel;
    elseif isa(channel, 'char') && strcmp(datasource, 'edf')
        edfchannel = channel;
    elseif length(channel) == 2 && strcmp(datasource, 'both')
        eegchannel = channel{1};
        edfchannel = channel{2};
    else
        error('Invalid channel and datasource inputs. Please check.')
    end
end
assert(isa(edfchannel, 'char'), 'Specified EDF channel must be a string for SleepEEG_loadedf() to work.')

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';

%%
disp('--------------------------')
disp([mHead 'SleepEEG_hypno()']);
disp('--------------------------')

tic

%%
% Add path to the necessary codes
SleepEEG_addpath(matlabroot, project);

% Construct correct paths and filenames
[dataDir, ~, ~, ~] = SleepEEG_configDir(subID, fnsuffix, false, project);

clinicalpath = fullfile(dataDir, subID, 'clinical');
analysispath = fullfile(dataDir, subID, 'analysis');
edfFN = [subID, '_', fnsuffix, '_clinical_deidentified.edf'];
if contains(fnsuffix, 'night1')
    txtFN = [subID, '_night1_scoring.txt'];
    alignedfFN = [subID, '_night1_aligned.edf'];
    rvsFN = [subID, '_night1_rvsalign.mat'];
elseif contains(fnsuffix, 'night2')
    txtFN = [subID, '_night2_scoring.txt'];
    alignedfFN = [subID, '_night2_aligned.edf'];
    rvsFN = [subID, '_night2_rvsalign.mat'];
else
    error('No night string is found in the fnsuffix. Please double check.')
end

tempstr = strsplit(txtFN, '_scoring.txt');
plotID = tempstr{1};

%% Extract the staging channel
[ stage_channel, EEGvec_stage_channel ] = SleepEEG_readstaging(txtFN, edfFN, alignedfFN ,clinicalpath, rvsFN);

%% Plot spectrograms
% two possible datasources to grab electrode data from
disp([mHead, 'Plotting spectrograms with the hypnogram overlayed...'])

if strcmp(datasource, 'eeg') ||  strcmp(datasource, 'both')
    %% Compute spectrogram based on an EEG channel from the EEG data
    % with hypnogram aligned and save the plot to "analysis"

    EEG = SleepEEG_loadset(subID, fnsuffix, project, true, 'Z3', true);
    % confirm that the length of EEGvec_stage_channel is correct
    assert(length(EEGvec_stage_channel) == size(EEG.data,2), 'Length of stage channel in EEG timeframe is incorrect.')

    % figure out which channel to compute the spectrogram on
    eegchannel_num = 0;
    if isa(eegchannel, 'char')
        for ii = 1:length(EEG.chanlocs)
            if strcmp(EEG.chanlocs(ii).labels, eegchannel)
                eegchannel_num = ii;
            end
        end
    else
        eegchannel_num = eegchannel;
    end
    assert(eegchannel_num > 0, 'Invalid channel number for the EEG .set data.')

    Fs = EEG.srate;
    t = EEG.times/1000;

    % calculate the spectrogram
    [mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_mex(EEG.data(eegchannel_num,:),...
        Fs, [0, Fs/2], [15, 29], [30, 5], 2^14, 'linear', 'unity', false, false);

    % plot spectrogram
    visfreq = [0,40];
    freq_idx = sfreqs>=visfreq(1) & sfreqs <= visfreq(2);

    f1 = figure;
    ax = figdesign(3,1,'merge',{2:3});
    for ii = 1:length(ax); title(ax(ii), ii); end
    set(f1, 'Units', 'inches');
    set(f1, 'Position', [1 1 16 10]);

    axes(ax(1))
    hypnoplot(t/60, EEGvec_stage_channel);
    title('Hypnogram from Sleep Scoring')
    set(gca, 'FontSize', 16)
    colorbar

    axes(ax(2))
    imagesc(stimes/60, sfreqs(freq_idx), pow2db(mt_spectrogram(freq_idx,:)))
    axis xy
    title(['Overnight Sleep - ', EEG.chanlocs(eegchannel_num).labels])
    xlabel('Times (min)')
    ylabel('Frequency (Hz)')
    colormap jet
    climscale;
    axis tight
    c = colorbar;
    ylabel(c, 'PSD (dB)');
    set(gca, 'FontSize', 16)

    linkaxes([ax(1), ax(2)], 'x')
    saveas(gcf, fullfile(analysispath, [plotID '_spectrohypnogram_EEG_',EEG.chanlocs(eegchannel_num).labels,'.png']))
    savefig(fullfile(analysispath, [plotID '_spectrohypnogram_EEG_',EEG.chanlocs(eegchannel_num).labels,'.fig']))
end

if strcmp(datasource, 'edf') ||  strcmp(datasource, 'both')
    %% Compute spectrogram based on an EEG channel from the aligned EDF data
    % with hypnogram aligned and save the plot to "clinical"

    % figure out which analog channel
    if strcmp(edfchannel, 'F3')
        analogchan = 'LL2';
    elseif strcmp(edfchannel, 'F4')
        analogchan = 'RR2';
    elseif strcmp(edfchannel, 'C3')
        analogchan = 'LA2';
    elseif strcmp(edfchannel, 'C4')
        analogchan = 'RA2';
    elseif strcmp(edfchannel, 'O1')
        analogchan = 'LL11';
    elseif strcmp(edfchannel, 'O2')
        analogchan = 'RR11';
    else
        error('Specified edfchannel is not valid. Please check.')
    end

    % grab a channel of signal in the aligned EDF+ file
    [header, signalHeader, signalCell] = SleepEEG_loadedf(fullfile(clinicalpath, alignedfFN), {edfchannel});
    aligned_signal = signalCell{1};
    % create a time vector
    Fs = signalHeader.samples_in_record / header.data_record_duration;
    t = 0:1/Fs:1/Fs*(length(aligned_signal)-1);

    % calculate the spectrogram
    [mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_mex(aligned_signal,...
        Fs, [0, Fs/2], [15, 29], [30, 5], 2^14, 'linear', 'unity', false, false);

    % plot spectrogram
    visfreq = [0,40];
    freq_idx = sfreqs>=visfreq(1) & sfreqs <= visfreq(2);

    f1 = figure;
    ax = figdesign(3,1,'merge',{2:3});
    for ii = 1:length(ax); title(ax(ii), ii); end
    set(f1, 'Units', 'inches');
    set(f1, 'Position', [1 1 16 10]);

    axes(ax(1))
    hypnoplot(t/60, stage_channel);
    title('Hypnogram from Sleep Scoring')
    set(gca, 'FontSize', 16)
    c=colorbar;
    set(c, 'yticklabel', [])

    axes(ax(2))
    imagesc(stimes/60, sfreqs(freq_idx), pow2db(mt_spectrogram(freq_idx,:)))
    axis xy
    title(['Overnight Sleep - ', edfchannel, ' analog (',analogchan, ')' ])
    xlabel('Times (min)')
    ylabel('Frequency (Hz)')
    colormap jet
    climscale;
    axis tight
    c = colorbar;
    caxis([-20, 20])
    ylabel(c, 'PSD (dB)');
    set(gca, 'FontSize', 16)

    linkaxes([ax(1), ax(2)], 'x')
    saveas(gcf, fullfile(clinicalpath, [plotID '_spectrohypnogram_EDF_',edfchannel,'.png']))
    savefig(fullfile(clinicalpath, [plotID '_spectrohypnogram_EDF_',edfchannel,'.fig']))
end

%%
close all
disp([mHead, 'SleepEEG_hypno() completed. Total time taken: '])
toc

end
