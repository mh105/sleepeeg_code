function [ data ] = inspect_channel(EEG1,channel, duration)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
data = EEG1.data(label2num(channel, EEG1.chanlocs), EEG1.times > (EEG1.times(end)-5*60*1000) &...
    EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
assert(length(data)/EEG1.srate/60 == duration, 'Incorrect length of eyes closed.')

figure
hold on
plot(EEG1.times, EEG1.data(label2num(channel, EEG1.chanlocs), :))
plot(EEG1.times(EEG1.times > (EEG1.times(end)-5*60*1000) &...
    EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000)), data)


method = 'std';
hf_crit = 4;
hf_pass = 25;
bb_crit = 4;
bb_pass = .1;
smooth_duration = 2;
verbose = false;
histogram_plot = false;

%Detect artifacts
[artifacts, hf_artifacts, bb_artifacts, y_high, y_broad] = EEG_detect_time_domain_artifacts(data,...
    EEG1.srate, method, hf_crit, hf_pass, bb_crit, bb_pass, smooth_duration, verbose, histogram_plot);

data(artifacts) = nan;

figure; plot(artifacts, 'LineWidth', 5)


Fs = EEG1.srate;

[mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_optimized(data,...
    Fs, [0, (Fs/2)], [2,3], [5,1], 2^16, 'linear', false, true);
[mt_spectrogram1, stimes1, sfreqs1] = multitaper_spectrogram_optimized(data,...
    Fs, [0, (Fs/2)], [2,3], [10,1], 2^16, 'linear', false, true);
[mt_spectrogram2, stimes2, sfreqs2] = multitaper_spectrogram_optimized(data,...
    Fs, [0, (Fs/2)], [4,7], [10,1], 2^16, 'linear', false, true);

y = nanmean(mt_spectrogram, 1);
y1 = nanmean(mt_spectrogram1, 1);
y2 = nanmean(mt_spectrogram2, 1);

% What effects do the parameter choices have on the spectral
% estimation?
figure
set(gcf, 'Position', [1601 70 1600 1168])
hold on
plot(sfreqs, pow2db(y), 'LineWidth', 2)
plot(sfreqs1, pow2db(y1), 'LineWidth', 2)
plot(sfreqs2, pow2db(y2), 'LineWidth', 2)
xlim([0, 40])
legend('TW=2, N=5sec', 'TW=2, N=10sec', 'TW=4, N=10sec', 'Location', 'best')
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
title('Multitaper Spectral Estimation')
set(gca, 'FontSize', 20)


end

