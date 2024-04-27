%% Testing of the ANT EEG physical min max
% This is a test script used to understand the deal with physical min/max
% values in the ANT EEG signal. What min/max values should we use, and when
% exceeded what does that indicate about the signal quality we are using
% in the .set files? 
%
% This impacts how Amanda's trigger alignment code handles EEG signals and
% determines whether we need to saturate out-of-bound values. 
%
% Anyways, here we go. 
%
close all
clear all

% ANT system EDF specifications in Amanda's trigger alignment code
phy_min = -83886;
phy_max = 83886;
dig_min = -32768;
dig_max = 32767;

% channels in EDF:
edfchannels = {'F3', 'F4', 'C3', 'C4', 'O1', 'O2', 'M1', 'M2'};
channels = {'LL2', 'RR2', 'LA2', 'RA2', 'LL11', 'RR11', 'LD6', 'RD6'};

night = '2';
rootdir = '/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/EEG/data/JFIT_F_41';

%% use readcnt to load .cnt file and check for min/max
% This is to make sure the out of min/max values are not due to some
% post-processing we did after calling pop_loadeep_v4()
filepath = fullfile(rootdir, 'raw');
cntEEG = ANT_interface_readcnt(['JFIT_F_41_night', night, '_Sleep.cnt'], filepath, [false, 0], true);

%%
% generate some pictures from the intact cnt data (not downsampled)
for ii = 1:length(channels)
    
    channel = channels{ii};
    
    figure
    hold on
    t = cntEEG.times;
    x = cntEEG.data(label2num(channel, cntEEG.chanlocs), :);
    plot(t, x)
    plot(xlim, [phy_max, phy_max], 'g', 'LineWidth', 2)
    plot(xlim, [phy_min, phy_min], 'g', 'LineWidth', 2)
    title([edfchannels{ii}, ' - ', channels{ii}], 'FontSize', 20)
    exthresh = x>phy_max |  x<phy_min;
    scatter(t(exthresh), x(exthresh), 'r')
    
    picname = fullfile(rootdir,'min_max_investigation', ['night', night, '_cnt_fs1000_', channel, '.png']);
    saveas(gca, picname);
    
    close all
end

%%
% report the out-of-bound values in each channel
format long; clc
for ii = 1:length(channels)
    channel = channels{ii};
    x = cntEEG.data(label2num(channel, cntEEG.chanlocs), :);
    exthresh = x>phy_max |  x<phy_min;
    disp([edfchannels{ii}, ' - ', channels{ii}])
    disp(unique(x(exthresh)))
end

%% OK! I think I have found the culprit!!!
% So, -83886/83886 indeed seem to be physical min/max on the exported .cnt
% files when read in without downsampling. They exceed the physical min/max
% threshold as a result of downsampling! 
%
%
%

%% Read in the downsampled set file and check that min/max are indeed exceeded 
filepath = fullfile(rootdir,'set');
EEG1 = ANT_interface_loadset(['JFIT_F_41_night', night, '_Sleep_ds500_Z3.set'], filepath);

signals = zeros(length(channels), size(EEG1.data, 2));
for ii = 1:length(channels)
    signals(ii,:) = EEG1.data(label2num(channels{ii}, EEG1.chanlocs), :);
end

%% 
% plot the raw signal time trace 
for ii = 1:size(signals, 1)
    
    figure
    hold on
    t = EEG1.times;
    plot(t, signals(ii, :))
    plot(xlim, [phy_max, phy_max], 'g', 'LineWidth', 2)
    plot(xlim, [phy_min, phy_min], 'g', 'LineWidth', 2)
    title([edfchannels{ii}, ' - ', channels{ii}], 'FontSize', 20)
    exthresh = signals(ii, :)>phy_max |  signals(ii, :)<phy_min;
    scatter(t(exthresh), signals(ii, exthresh), 'r')
    
    picname = fullfile(rootdir,'min_max_investigation', ['night', night, '_set_fs500_', channels{ii}, '.png']);
    saveas(gca, picname);
    
    % detrend
    figure
    hold on
    t = EEG1.times;
    tempsignal = detrend(signals(ii, :));
    plot(t, tempsignal)
    plot(xlim, [phy_max, phy_max], 'g', 'LineWidth', 2)
    plot(xlim, [phy_min, phy_min], 'g', 'LineWidth', 2)
    title(['Detrended ', edfchannels{ii}, ' - ', channels{ii}], 'FontSize', 20)
    exthresh = tempsignal>phy_max |  tempsignal<phy_min;
    scatter(t(exthresh), tempsignal(exthresh), 'r')
    
    picname = fullfile(rootdir,'min_max_investigation', ['night', night, '_set_fs500_', channels{ii}, '_detrend.png']);
    saveas(gca, picname);
    
    close all
end

%%
% Ok, so set file is essentially the same as the downsampled .cnt read
% using the ANT_interface_readcnt() function. Fine. So the thresholds are
% exceeded as a result of downsampling. 
%
% 
%
% Now the tricky question is: what am I seeing in the exported EDF files?
% What the hell are these signals??? Did I just export things wrong?

%% Compare against ANT exported EDF
edfname = fullfile(rootdir, 'raw/', 'JFIT_F__41_2019-11-11_22-20-32_Segment_0.edf');
[header_all, signalHeader_all] = SleepEEG_loadedf(edfname);

channel_to_load = {};
for ii = 1:length(channels)
   channel_to_load{ii} = ['EEG ' channels{ii} '-Z3']; 
end

% load the signals 
[header, signalHeader, signalCell] = SleepEEG_loadedf(edfname, channel_to_load);
edf_FS = double(signalHeader(1).samples_in_record) / double(header.data_record_duration);

for ii = 1:length(signalCell)
    edf_channel = signalCell{ii};
    edf_time = linspace(1, length(edf_channel)/edf_FS, length(edf_channel));
    
    figure
    plot(edf_time, edf_channel)
    title([edfchannels{ii}, ' - ', channels{ii}], 'FontSize', 20)
    
    picname = fullfile(rootdir,'min_max_investigation', ['night', '1', '_edf_fs1000_', channels{ii}, '.png']);
    saveas(gca, picname);
    
    close all
end

% something is wrong about these exported .edf files. We need to re-export
% and check! 

%% What's the actual impact of using hard-coded -83886/83886 min/max? 

EEG1 = SleepEEG_loadset('TDelp', 'resting');

%%
close all
channel = 25;
deamp_factor = 20;
x = EEG1.data(channel,:)/deamp_factor;

figure
hold on
t = EEG1.times;
plot(t, x)
plot(xlim, [phy_max, phy_max], 'g', 'LineWidth', 2)
plot(xlim, [phy_min, phy_min], 'g', 'LineWidth', 2)
title(EEG1.chanlocs(channel).labels, 'FontSize', 20)
exthresh = x>phy_max |  x<phy_min;
scatter(t(exthresh), x(exthresh), 'r')

% spectral estimates 
% calculate the spectrogram
[mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_optimized(x, EEG1.srate, [0, (EEG1.srate/2)]);

% plot spectrogram
visfreq = [0,70];
freq_idx = sfreqs>=visfreq(1) & sfreqs <= visfreq(2);

figure
imagesc(stimes, sfreqs(freq_idx), pow2db(mt_spectrogram(:, freq_idx)'))
axis xy
title('Resting Eyes Closed - Z10')
xlabel('Times (s)')
ylabel('Frequency (Hz)')
colormap jet
climscale
axis tight
c= colorbar;
ylabel(c, 'PSD (dB)');
set(gca, 'FontSize', 16)

% plot the spectrum
y = mean(mt_spectrogram, 1);
h = figure;
hold on
plot(sfreqs, pow2db(y), 'LineWidth', 2)
title('Mean Spectrum - Z10')
xlabel('Frequency (Hz)')
ylabel('Average PSD (dB)')
set(gca, 'FontSize', 16)
axis tight
xlim(visfreq)
ylim(ylim*1.1)

% Now we convert this channel using physical min/max to digital min/max and
% write with 16bit integers and back.

signal = (x-phy_min)/(phy_max-phy_min);
signal = signal.*double(dig_max-dig_min)+dig_min;

figure
hold on
t = EEG1.times;
plot(t, signal)

signal_int16 = int16(signal);
plot(t, signal_int16)

newx = double(signal_int16);
% Convert from digital to physical values
value = (newx-dig_min)/(dig_max-dig_min);
value = value.*double(phy_max-phy_min)+phy_min;

figure
hold on
plot(t,x)
plot(t,value)

% spectral estimate
[mt_spectrogram1, stimes1, sfreqs1] = multitaper_spectrogram_optimized(value, EEG1.srate, [0, (EEG1.srate/2)]);

% plot spectrogram
visfreq = [0,70];
freq_idx = sfreqs1>=visfreq(1) & sfreqs1 <= visfreq(2);

figure
imagesc(stimes1, sfreqs1(freq_idx), pow2db(mt_spectrogram1(:, freq_idx)'))
axis xy
title('Resting Eyes Closed - Z10')
xlabel('Times (s)')
ylabel('Frequency (Hz)')
colormap jet
climscale
axis tight
c= colorbar;
ylabel(c, 'PSD (dB)');
set(gca, 'FontSize', 16)

% plot the spectrum
y1 = mean(mt_spectrogram1, 1);
figure(h)
hold on
plot(sfreqs, pow2db(y1), 'LineWidth', 2)
title(['Mean Spectrum - Z10 Comparison to shrinkage factor ', num2str(deamp_factor)])
xlabel('Frequency (Hz)')
ylabel('Average PSD (dB)')
legend('Original', 'Converted')
set(gca, 'FontSize', 16)
axis tight
xlim(visfreq)
ylim(ylim*1.1)











