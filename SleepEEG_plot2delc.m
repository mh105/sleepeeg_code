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

%% Configure inputs 

assert(exist('chanlocs', 'var') || exist('EEG', 'var'), 'Either chanlocs or EEG structure needs to be inputed for now.')

if ~exist('EEG', 'var')
    EEG = eeg_emptyset();
    EEG.chanlocs = chanlocs;    
end

%%
% visualize electrode locations 
f = figure;
set(f, 'Units', 'inches');
set(f, 'Position', [1 1 16 10]);

topoplot([],EEG.chanlocs,'style','both','electrodes','ptsnumbers','emarker', {'.', 'k', 15, 1});

L = findobj(gcf, 'type', 'Text');
for ind = 1:length(EEG.chanlocs)
    set(L(length(EEG.chanlocs)+1-ind), 'FontSize', 14)
end

title([subID ' Channel Locations'], 'FontSize', 30, 'Interpreter', 'none')

end
