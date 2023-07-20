function [ ] = anabehav_all(subID)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - ALL**
%
% - Preprocess all cognitive task data in the AD Sleep project: 
%   1. Motor Sequence Task (MST)
%   2. Psychomotor Vigilance Task (PVT)
%   3. Auditory Oddball Task (P300)
%   4. Object Location Task (ObjLoc)
%   5. Process Dissociation Paradigm (PDP)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for
%                           relevant task txt files in the 'task' folder
%                           Filenames are hardcoded in the analysis
%                           functions.  
%
%                           if any of required files is missing, an error
%                           will be reported and the function will
%                           terminate.  
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Addpath to appropriate folders
% matlabroot will help detect the current environment
dataDir = SleepEEG_addpath(matlabroot);

%% Iterate through all tasks 
filepath = fullfile(dataDir, subID, 'task');

anabehav_MST(subID, filepath);
anabehav_PVT(subID, filepath);
anabehav_P300(subID, filepath);
anabehav_OBJLOC(subID, filepath);
anabehav_PDP(subID, filepath);

end




