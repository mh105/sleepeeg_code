function [ EB ] = SleepEEG_bridge(subID, fnsuffix, project, EEG, fileID, outputDir)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - BRIDGE**
%
% - checks bridged electrodes during the EEG recording
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
%           - fileID:       name of the .cnt file (no .cnt suffix).
%
%           - outputDir:    directory path to folder for dumping all
%                           analyses outputs.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - EB:           a structure containing various bridging
%                           analysis outputs. This is also saved into the
%                           analysis folder with the _bridging_EB_structure
%                           suffix.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';

%%
disp('---------------------------')
disp([mHead 'SleepEEG_bridge()']);
disp('---------------------------')

%% Load Data if the EEG structure is not an input
if nargin < 4
    [ dataDir, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
    assert(isfile(fullfile(dataDir, subID, 'set', [subID, '_', fnsuffix, '_ds500_Z3.set'])),...
        [mHead, 'Downsampled data set .set file is not available for ' fileID])
    EEG = SleepEEG_loadset(subID, fnsuffix, project);
end

if ~exist('fileID', 'var')
    [ ~, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

% Bridging detection isn't affected by re-referencing, so we won't append
% this information to the saved file name.
% fileID = [fileID, '_', EEG.refscheme];

%% Bridging electrode detection
% Use eBridge based on temporal variance on pairwise comparisons between
% channels to check bridigng electrodes. Only applied between the first and
% last occurencess of time-locking events padded by the epoch limits - to
% reduce the processing time of eBridge.m
tic

% Calls eBridge.m
EB = eBridge(EEG, 'PlotMode', 2, 'BCT', 0.5, 'EpLength', EEG.srate, 'FiltMode', 1);
set(gca, 'FontSize', 16)
title([fileID ' ED distribution'], 'Interpreter', 'none')
saveas(gcf, fullfile(outputDir, [fileID '_bridging_analysis.png']))

% Optional: visualize bridged electrodes
figure; topoplot([],EEG.chanlocs,'style','both','electrodes','ptslabels','emarker', {'.', 'k', 15, 1});
L = findobj(gcf, 'type', 'Text');
for ind = EB.Bridged.Indices
    set(L(length(EEG.chanlocs)+1-ind), 'FontSize', 20)
    set(L(length(EEG.chanlocs)+1-ind), 'Color', [1 0 0])
end
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1:2) pos(3:4)*1.5])
title([fileID ' Bridged'], 'FontSize', 30, 'Interpreter', 'none')
saveas(gcf, fullfile(outputDir, [fileID '_bridged_electrodes.png']))

close all

%% Save the bridging analysis output file
save(fullfile(outputDir, [fileID '_bridging_EB_structure']), 'EB', '-v7.3')

disp([mHead, 'Total time taken in bridging analysis...']);
disp(' ')
toc

end
