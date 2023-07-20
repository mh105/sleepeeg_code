function [] = SleepEEG_plot2delc(subID, chanlocs, EEG)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - PLOT2DELC**
%
% - plots the 2D topoplot of EEG electrodes 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - chanlocs:     a EEGLAB chanloc structure 
%
%           - EEG:          data structure containing the EEG data and
%                           chanloc.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%% Configure inputs 

assert(exist('chanlocs', 'var') || exist('EEG', 'var'), 'Either chanlocs or EEG structure needs to be inputted for now.')

if ~exist('EEG', 'var')
    EEG = eeg_emptyset();
    EEG.chanlocs = chanlocs;    
end

%%
% visualize electrode locations 
figure
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1:2) pos(3:4)*1.5])

topoplot([],EEG.chanlocs,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});

L = findobj(gcf, 'type', 'Text');
for ind = 1:length(EEG.chanlocs)
    set(L(length(EEG.chanlocs)+1-ind), 'FontSize', 14)
end

title([subID ' Channel Locations'], 'FontSize', 30, 'Interpreter', 'none')

end
