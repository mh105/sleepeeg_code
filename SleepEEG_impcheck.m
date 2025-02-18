function [ High_imp_electrode ] = SleepEEG_impcheck(subID, fnsuffix, project, EEG, fileID, outputDir, visualize)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - IMPCHECK**
%
% - checks impedances at the beginning and the end of EEG recording
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
%           - visualize:    whether to generate plots of impedances
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - High_imp_electrode:
%                           a vector of double type containing the channel
%                           number of the electrodes with high impedance
%                           values.
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
disp('-----------------------------')
disp([mHead 'SleepEEG_impcheck()']);
disp('-----------------------------')

%% Load Data if the EEG structure is not an input
if nargin < 4
    [ dataDir, ~, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, false, project);
    assert(isfile(fullfile(dataDir, subID, 'set', [subID, '_', fnsuffix, '_ds500_Z3.set'])),...
        [mHead, 'Downsampled data set .set file is not available for ' fileID])
    EEG = SleepEEG_loadset(subID, fnsuffix, project);
end

if ~exist('visualize', 'var')
    visualize = true;
end

%% Plot impedance distributions
% mask dead electrodes
if isempty(EEG.endimp) % for some early asalab recordings, only initial impedance, no end impedance measures
    EEG.endimp = EEG.initimp;
end
EEG.initimp(EEG.initimp > 500) = NaN;
EEG.endimp(EEG.endimp > 500) = NaN;

if visualize
    % Save a .fig of the electrode location topoplot
    SleepEEG_plot2delc(subID, EEG.chanlocs)
    savefig(fullfile(outputDir, [fileID '_montage_topoplot.fig']))
    
    % Initial Impedance
    figure
    histogram(EEG.initimp, 50)
    xlabel('Impedance k\Omega')
    ylabel('count')
    title([fileID ' Initial Imps'], 'Interpreter', 'none')
    set(gca, 'FontSize', 16)
    saveas(gcf, fullfile(outputDir, [fileID '_InitImp_dist.png']))
    
    % End Impedance
    figure
    histogram(EEG.endimp, 50)
    xlabel('Impedance k\Omega')
    ylabel('count')
    title([fileID ' End Imps'], 'Interpreter', 'none')
    set(gca, 'FontSize', 16)
    saveas(gcf, fullfile(outputDir, [fileID '_EndImp_dist.png']))
    
    close all
end

%% Report electrodes with impedance above 50kOhm by the end of recording
High_imp_electrode = [find(EEG.endimp > 50); EEG.endimp(EEG.endimp > 50)];

if ~isempty(High_imp_electrode)
    disp([mHead, 'High impedance electrodes at the end of recording:'])
    fprintf([mSpace, 'Channel Number:'])
    disp(High_imp_electrode(1,:))
    fprintf([mSpace, 'Impedance value:'])
    disp(High_imp_electrode(2,:))
    disp(' ')
else
    disp('No high impedance electrode detected.')
end

end
