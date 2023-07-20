function [] = quick_spectrogram(mt_spectrogram, stimes, sfreqs)

% plot spectrogram
visfreq = [0,40];
freq_idx = sfreqs >= visfreq(1) & sfreqs <= visfreq(2);
figure
imagesc(stimes, sfreqs(freq_idx), nanpow2db(mt_spectrogram(:, freq_idx)'))
axis xy
title('Spectrogram')
xlabel('Times (s)')
ylabel('Frequency (Hz)')
colormap jet
climscale
axis tight
c= colorbar;
ylabel(c, 'Power (dB)');
set(gca, 'FontSize', 20)

end