function [ Deadidx ] = SleepEEG_dead(subID, fnsuffix, project, EEG, fileID, outputDir, visualize)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - DEAD**
%
% - Used to detect dead electrodes and electrodes with high impedances
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
%           - visualize:    whether to plot a topoplot of dead electrodes
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - Deadidx:      a vector of double type containing the channel
%                           number of electrodes detected to be dead. The
%                           reference electrode is excluded from this list
%                           because it has constant zero recording by
%                           definition.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('-------------------------')
disp([mHead 'SleepEEG_dead()']);
disp('-------------------------')

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

if ~exist('visualize', 'var')
    visualize = true;
end

%% Dead electrode detection
% flag dead electrodes by an impedance threshold
if isempty(EEG.endimp) % for some early asalab recordings, only initial impedance, no end impedance measures
    EEG.endimp = EEG.initimp;
end
EEG.initimp(EEG.initimp > 500) = NaN;
endimp_copy = EEG.endimp;
EEG.endimp(EEG.endimp > 500) = NaN;

% Visualize dead electrodes as long as we have end impedance measures
if ~isempty(EEG.endimp)

    Deadidx = [];
    for i = 1:size(EEG.data, 1)
        % Flagged dead either initial or end impedance > 500kOhm
        if (~isempty(EEG.initimp) && isnan(EEG.initimp(i))) || isnan(EEG.endimp(i))
            % crude detection at the moment, a better way will be to segment
            % the data and compute whether it is dead for each segment
            % separately to catch changes over the night
            Deadidx(length(Deadidx)+1) = i; %#ok<AGROW>
            %         elseif nanmean(EEG.data(i,:)) == 0 % asalab dead channels have 0 mV
            %             Deadidx(length(Deadidx)+1) = i;
        end
    end

    if visualize
        % visualize dead electrodes
        figure; topoplot([],EEG.chanlocs,'style','both','electrodes','ptslabels','emarker', {'.', 'k', 15, 1});
        L = findobj(gcf, 'type', 'Text');
        for ind = Deadidx
            set(L(length(EEG.chanlocs)+1-ind), 'FontSize', 20)
            set(L(length(EEG.chanlocs)+1-ind), 'Color', [0 1 0])
        end
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1:2) pos(3:4)*1.5])
        if ~isempty(EEG.initimp)
            title_text = [fileID ' Dead (Init & End)'];
        else
            title_text = [fileID ' Dead (End Only)'];
        end
        title(title_text, 'FontSize', 30, 'Interpreter', 'none')
        saveas(gcf, fullfile(outputDir, [fileID '_dead_electrodes.png']))

        close all
    end

else
    warning('End Impedance values are unavailable!')
    warning('This might be because the recording was not ended properly, likely from software crashed half way!')
    Deadidx = [];

end

%% Remove the unipolar reference electrode from outputted Deadidx
ref_channel = nan;
for ii = 1:length(EEG.chanlocs)
    if strcmp(EEG.refscheme, EEG.chanlocs(ii).labels)
        ref_channel = ii;
    end
end
Deadidx(Deadidx == ref_channel) = [];

%% Report dead electrodes
if ~isempty(Deadidx)
    disp([mHead, 'Dead electrodes in the recording:'])
    fprintf([mSpace, 'Channel Number:'])
    disp(Deadidx)
    fprintf([mSpace, 'Impedance value:'])
    disp(endimp_copy(Deadidx))
    disp(' ')
else
    disp('No dead electrodes detected.')
end

end
