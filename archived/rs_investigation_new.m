%% Resting State Spectral Estimation Investigation
% Investigate the Resting State Data in the 5 subjects we have:
% - What referencing should we use?
% - What detrending should we use?
% - What spectral estimation method should we use? 
% - What multitaper parameters should we use? 
% - Confirm that the spectral estimations are consistent within the same
% session between Z2 and Z10 channels for the high frequency noise and low
% frequnecy slow oscillation but different on alpha band. 
% - Confirm that the spectral estimations are consistent within a subject
% across the two nights but different across subjects

close all
clear all
clc

% Change the current folder to the folder of this m-file.
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))
clearvars tmp

dataDir = SleepEEG_addpath(matlabroot);

% addpath to artifact detection folder
addpath(genpath('/Users/alexhe/Dropbox_Personal/Mac_Desktop/MIT_HST/NSRL_Purdon_Prerau/artifact_detection'))

%% Load all subjects' data
subject_list = {'JALL_M_00', 'KROA_F_90', 'TTRA_M_89', 'JFIT_F_41', 'MDEM_F_42'};
suffix_list = {'night1_Resting', 'night2_Resting'};

mean_alpha = table;
mean_alpha.subID = string;
mean_alpha.task = string; 
EEG_cell = {};
EEG_apeak = {};
Fs = 500;

for ii = 1:length(subject_list)
    
    sid = subject_list{ii};
    mean_alpha.subID((ii-1)*2+1) = sid;
    mean_alpha.subID((ii-1)*2+2) = sid;
    
    for j = 1:length(suffix_list)
        
        duration = 4.5;
        
        task = suffix_list{j};
        mean_alpha.subID((ii-1)*2+j) = task;
        
        % load the set file
        EEG1 = SleepEEG_loadset(sid, task);
        
        % get the Z10 channel referenced to Z3
        data_z3 = EEG1.data(label2num('Z10', EEG1.chanlocs), EEG1.times > (EEG1.times(end)-5*60*1000) &...
            EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
        
        % exclude dead electrodes and the EOG electrode
        deadidx = find((EEG1.endimp > 500 | EEG1.initimp > 500) == 1);
        EEG1.data([deadidx, 129], :) = [];
        EEG1.chanlocs([deadidx, 129]) = [];
        
        % get the common average channel using remaining channels
        ca_channel = mean(EEG1.data, 1);
        data_ca = ca_channel(EEG1.times > (EEG1.times(end)-5*60*1000) &...
            EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
        
        % re-reference to common average
        EEG1.data = EEG1.data - ca_channel;
        % store the EEG structure
        EEG_cell{(ii-1)*2+j,1} = EEG1;
        
        % grab eyes closed data from channel Z10 from last 4.5min of recording, excluding
        % the last 30 seconds due to movement artifacts at the end
        data = EEG1.data(label2num('Z10', EEG1.chanlocs), EEG1.times > (EEG1.times(end)-5*60*1000) &...
            EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
        assert(length(data)/Fs/60 == duration, 'Incorrect length of eyes closed.')
        
        
%         inspect_channel(EEG1, 'LE3', 2.5);

        
        % get channel Z2 as well for comparison
        data2 = EEG1.data(label2num('Z2', EEG1.chanlocs), EEG1.times > (EEG1.times(end)-5*60*1000) &...
            EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
        

        % visualize the time trace
        t = linspace(0,length(data)/Fs/60, length(data));
        
        figure
        set(gcf, 'Position', [1601 70 1600 1168])
        subplot(4,1,1)
        plot(t, data, 'LineWidth', 1)
        hold on
        plot(t, data_z3-data_ca, 'LineWidth', 1)
        title('Z10-CA Signal Time Trace')
        xlabel('Time (min)')
        ylabel('Voltage (\mu V)')
        set(gca, 'FontSize', 20)

        subplot(4,1,3)
        plot(t, data_z3, 'LineWidth', 1)
        title('Z10-Z3 Signal Time Trace')
        xlabel('Time (min)')
        ylabel('Voltage (\mu V)')
        set(gca, 'FontSize', 20)
        
        subplot(4,1,4)
        plot(t, data_ca, 'LineWidth', 1)
        title('CA-Z3 Signal Time Trace')
        xlabel('Time (min)')
        ylabel('Voltage (\mu V)')
        set(gca, 'FontSize', 20)
        
        % Artifact rejection on the data channel
        
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
            Fs, method, hf_crit, hf_pass, bb_crit, bb_pass, smooth_duration, verbose, histogram_plot);

        data(artifacts) = nan;
        
        subplot(4,1,2)
        plot(t, artifacts, 'LineWidth', 2)
        title('Z10-CA Artifact Detection')
        xlabel('Time (min)')
        ylabel('Artifact Boolean')
        set(gca, 'FontSize', 20)
        
        
        % save plot
        savefig(fullfile(dataDir, 'outputs', [sid, '_', task, '_timetrace']))
        saveas(gcf, fullfile(dataDir, 'outputs', [sid, '_', task, '_timetrace.png']))

%--------------------------------------------------------------------------     
%         % generate detrended data and visualize segment by segment
%         N = 10;
%
%         m = 11;
%
%         data_segment = data((m-1)*(N-1)*Fs+1:((m-1)*(N-1)+N)*Fs);
%         segment_t = t((m-1)*(N-1)*Fs+1:((m-1)*(N-1)+N)*Fs);
%         figure
%         subplot(2,1,1)
%         plot(segment_t, data_segment)
%         title(['Undetrended - ', num2str(N), 's'])
%         xlabel('Time (min)')
%         set(gca, 'FontSize', 20)
%         subplot(2,1,2)
%         plot(segment_t, detrend(data_segment))
%         title(['Detrended - ', num2str(N), 's'])
%         xlabel('Time (min)')
%         set(gca, 'FontSize', 20)

        %%
        % Investigate the different ways of using multitaper spectral
        % estimation on Z10-CA signal
        
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
        title({'Multitaper Spectral Estimation', [sid, '_', task]}, 'interpreter', 'none')
        set(gca, 'FontSize', 20)
        
        if j == 1
            ylimit = ylim;
        elseif j == 2
            ylim(ylimit)
        end
        
        % Ok, we will go with TW=2, N=5sec
        
        savefig(fullfile(dataDir, 'outputs', [sid, '_', task, '_meanSpectrum']))
        saveas(gcf, fullfile(dataDir, 'outputs', [sid, '_', task, '_meanSpectrum.png']))
        
        %% Topoplot of alpha peak power under CA
        chanlocs = EEG1.chanlocs;
        alpha_peakp = zeros(1, length(chanlocs));
        
        for k = 1:length(chanlocs)
            disp(['Processing channel #', num2str(k)])
            % grab channel
            s = EEG1.data(k, EEG1.times > (EEG1.times(end)-5*60*1000) &...
                EEG1.times <= (EEG1.times(end)-(5-duration)*60*1000));
            assert(length(s)/Fs/60 == duration, 'Incorrect length of eyes closed.')
            %Detect artifacts
            [artifacts, hf_artifacts, bb_artifacts, y_high, y_broad] = EEG_detect_time_domain_artifacts(s,...
                Fs, method, hf_crit, hf_pass, bb_crit, bb_pass, smooth_duration, verbose, histogram_plot);
            s(artifacts) = nan;
            % mtm one particular parameter set choice
            [mts, st, sf] = multitaper_spectrogram_optimized(s,...
                Fs, [0, (Fs/2)], [2,3], [5,1], 2^16, 'linear', false, false);
            % find alpha peak power
            y = nanmean(mts, 1);
            alpha_peakp(k) = max(y(sf >= 8 & sf <= 12));
        end
        
        EEG_apeak{ii, j, 1} = alpha_peakp;
        
        figure
        set(gcf, 'Position', [1601 70 1600 1168])
        topoplot(pow2db(alpha_peakp), chanlocs, 'style','both','electrodes','on','emarker', {'.', 'w', 30, 1});
        title({'8-12Hz Peak Power During Eyes Closed', 'Common Average Reference'}, 'FontSize', 16)
        c=colorbar;
        ylabel(c,'Power (dB)');
        set(gca, 'FontSize', 20)
        climscale
        
        if j == 1
            cspect_CA = caxis;
        elseif j == 2
            caxis(cspect_CA)
        end
        
        saveas(gcf, fullfile(dataDir, 'outputs', [sid, '_', task, '_topo_CA.png']))
        
        %% Topoplot of alpha peak power under Laplacian (Local Average)
        % load the set file
        EEG2 = SleepEEG_loadset(sid, task);
        
        % exclude the EOG electrode
        EEG2.data(129, :) = [];
        EEG2.chanlocs(129) = [];
        chanlocs = EEG2.chanlocs;
        
        % detect dead electrodes 
        deadidx = find((EEG2.endimp > 500 | EEG2.initimp > 500) == 1);
        
%         % load bridged electrodes
%         EB = load(fullfile(dataDir, sid, 'analysis', [sid, '_', task, '_bridging_EB_structure.mat']));
%         
        % get nearest neighbor Laplacian matrix
        load('duke_128_refmatrix.mat')
        L = squeeze(ref_matrix(end, :, :));
        
        % Loop through the matrix to modify local averages containing dead
        % electrodes 
        for k = 1:length(chanlocs)
            neighbors = find(L(k, :) ~= 0);
            % exclude the channel itself 
            neighbors(neighbors == k) = [];
            % check whether any of the neighbors is dead
            survivors = neighbors(~ismember(neighbors, deadidx));
            deadones = neighbors(ismember(neighbors, deadidx));
            if length(neighbors) ~= length(survivors)
                L(k, deadones) = 0;
                L(k, survivors) = -1 / length(survivors);
            end
        end
        
        % Laplacian re-reference 
        EEG2.data = L * EEG2.data;
        EEG2.refscheme = 'LP';
        
        % exclude dead electrodes
        EEG2.data(deadidx, :) = [];
        EEG2.chanlocs(deadidx) = [];
        chanlocs = EEG2.chanlocs;

        % generate topoplot
        alpha_peakp = zeros(1, length(chanlocs));
        
        for k = 1:length(chanlocs)
            disp(['Processing channel #', num2str(k)])
            % grab channel
            s = EEG2.data(k, EEG2.times > (EEG2.times(end)-5*60*1000) &...
                EEG2.times <= (EEG2.times(end)-(5-duration)*60*1000));
            assert(length(s)/Fs/60 == duration, 'Incorrect length of eyes closed.')
            %Detect artifacts
            [artifacts, hf_artifacts, bb_artifacts, y_high, y_broad] = EEG_detect_time_domain_artifacts(s,...
                Fs, method, hf_crit, hf_pass, bb_crit, bb_pass, smooth_duration, verbose, histogram_plot);
            s(artifacts) = nan;
            % mtm one particular parameter set choice
            [mts, st, sf] = multitaper_spectrogram_optimized(s,...
                Fs, [0, (Fs/2)], [2,3], [5,1], 2^16, 'linear', false, false);
            % find alpha peak power
            y = nanmean(mts, 1);
            alpha_peakp(k) = max(y(sf >= 8 & sf <= 12));
        end
        
        EEG_apeak{ii, j, 2} = alpha_peakp;
        
        figure
        set(gcf, 'Position', [1601 70 1600 1168])
        topoplot(pow2db(alpha_peakp), chanlocs, 'style','both','electrodes','on','emarker', {'.', 'w', 30, 1});
        title({'8-12Hz Peak Power During Eyes Closed', 'Laplacian Reference'}, 'FontSize', 16)
        c=colorbar;
        ylabel(c,'Power (dB)');
        set(gca, 'FontSize', 20)
        climscale
        
        if j == 1
            cspect_LP = caxis;
        elseif j == 2
            caxis(cspect_LP)
        end
        
        saveas(gcf, fullfile(dataDir, 'outputs', [sid, '_', task, '_topo_LP.png']))
        
    end
end

%%


























