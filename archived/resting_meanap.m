function [] = resting_meanap(subID,fnsuffix, interpolate)
%Function used to compute the mean alpha psd
if nargin < 3
    interpolate = true;
end

[dataDir, ~, fileID, ~] = SleepEEG_configDir(subID, fnsuffix, true);

% check whether the common averaged data already exists
savefn = [subID, '_', fnsuffix, '_ds500_AR.set'];

% First, we need to re-reference the data to common average reference
if isfile(fullfile(dataDir, subID, 'set', savefn))
    EEG = ANT_interface_loadset(savefn, fullfile(dataDir, subID, 'set'), true, true);
else
    % First, we need to re-reference the data to common average reference
    EEG = SleepEEG_call(subID, fnsuffix, 'none', 'downsample', [true, 500],...
        'dead', true, 'interp', interpolate, 'reref', {true, 'AR'}, 'oversave', true);
end

% Now we have the data, we compute the mean alpha psd and plot
% spectrogram - pick channel 25 to compute.

data = EEG.data(25, EEG.times > (EEG.times(end)-5*60*1000));

% calculate the spectrogram
[mt_spectrogram, stimes, sfreqs] = multitaper_spectrogram_optimized(data, EEG.srate, [0, (EEG.srate/2)]);
close all

% plot spectrogram
visfreq = [0,40];
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
saveas(gcf, fullfile(dataDir, subID, 'analysis', [fileID, '_meanAlpha_spectrogram.png']))

% plot the spectrum
y = mean(mt_spectrogram, 1);
figure
plot(sfreqs, pow2db(y), 'LineWidth', 2)
title('Mean Spectrum - Z10')
xlabel('Frequency (Hz)')
ylabel('Average PSD (dB)')
set(gca, 'FontSize', 16)
axis tight
xlim(visfreq)
ylim(ylim*1.1)
saveas(gcf, fullfile(dataDir, subID, 'analysis', [fileID, '_meanAlpha_spectrum.png']))

% compute the mean alpha PSD between 8-12
meanap = pow2db(mean(y(sfreqs>=8 & sfreqs <=12)));
disp(['Mean alpha psd = ', num2str(meanap) ,' dB'])
f = fopen(fullfile(dataDir, subID, 'analysis', [fileID, '_meanAlpha_value.txt']), 'w');
fprintf(f, '%.5f', meanap);
fclose(f);

% Done!
end

