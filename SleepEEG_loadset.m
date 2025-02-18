function [ EEG ] = SleepEEG_loadset(subID, fnsuffix, project, downsampled, reref, verbose)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - LOADSET**
%
% - Used to load .set EEG files based on subject ID and file suffix
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep".
%
%                           ***the final file name is in the form: subID_fnsuffix_ds500.set***
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
%           - downsampled:  whether to add '_ds500' when loading data.
%                           default: true
%
%           - reref:        a string specifying the reference scheme of the
%                           data to be loaded.
%                           default: 'Z3'
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - EEG:          an EEGLAB structure containing all information
%                           of the recording in .set file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    project = '';
    downsampled = true;
    reref = 'Z3';
    verbose = true;
elseif nargin < 4
    downsampled = true;
    reref = 'Z3';
    verbose = true;
elseif nargin < 5
    reref = 'Z3';
    verbose = true;
elseif nargin < 6
    verbose = true;
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';

%%
disp('----------------------------')
disp([mHead 'SleepEEG_loadset()']);
disp('----------------------------')

%% Define Directories of Codes and Data Folders
[dataDir, ~, ~, ~] = SleepEEG_configDir(subID, fnsuffix, verbose, project);

%% Configure filename and filepath
% for now, filename is hard-coded to load downsampled .set file, we can
% make it an argument input to load different ones.
if downsampled
    filename = [subID, '_', fnsuffix, '_ds500_', reref, '.set'];
else
    filename = [subID, '_', fnsuffix, '_', reref, '.set'];
end
filepath = fullfile(dataDir, subID, 'set');

%% Load .set file
EEG = ANT_interface_loadset(filename, filepath, verbose, false);

end
