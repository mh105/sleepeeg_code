function [ dataDir, datafn, fileID, outputDir ] = SleepEEG_configDir(subID, fnsuffix, verbose)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - CONFIGDIR**
%
% - used to configure directories and specify path variables. Since lots of
% SleepEEG functions require to specify the right directories and places to
% find data, these actions are consolidated in this function
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep". if specified as
%                           [] empty, then all .cnt files under the subID
%                           directory will be returned as a cell array for
%                           "fileID" output.
%
%                           ***the final file name is in the form: subID_fnsuffix.cnt***
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Outputs:
%
%           - dataDir:      directory path containing all subjects' data.
%
%           - datafn:       full path to the .cnt file including preceding
%                           directories.
%
%           - fileID:       name of the .cnt file (no .cnt suffix). if
%                           fnsuffix input is specified as [] empty, then
%                           all .cnt files under the subID directory will
%                           be returned in datafn and subjectId as a cell
%                           array.
%
%           - outputDir:    directory path to folder for dumping all
%                           analyses outputs.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    verbose = true;
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
if verbose
    disp('------------------------------')
    disp([mHead 'SleepEEG_configDir()']);
    disp('------------------------------')
end

%% Addpath to appropriate folders
% matlabroot will help detect the current environment
dataDir = SleepEEG_addpath(matlabroot);

%% Instantiate Subject Data and Path Directories

if isempty(fnsuffix) % if fnsuffix is not specified, we process all available data for the subject
    cntfiles = dir(fullfile(dataDir, subID, 'raw', [subID, '*.cnt']));
    if verbose
        disp("===================================================")
        disp([mHead, 'Only subject ID is provided to SleepEEG_configDir(), ALL available .cnt files will be processed!'])
        disp([mHead, '.cnt files under "', subID, '" folder to be processed are:'])
    end
    datafn = []; % returns an empty data filename since no single file
    fileID = cell(1,length(cntfiles));
    % Now 'fileId' is a 1 x m vector with fnsuffix for each .cnt files
    
    for i = 1:length(cntfiles)
        C = strsplit(cntfiles(i).name, {subID, '.cnt'});
        if C{2}(1) == '_'
            fileID{i} = C{2}(2:end);%#ok<*AGROW>
        else
            fileID{i} = C{2};
        end
        if verbose; disp(cntfiles(i).name); end
    end
    if verbose; disp(" "); end
    
else % otherwise use fnsuffix to select one .cnt file
    cntfn = [subID, '_', fnsuffix, '.cnt'];
    datafn =  fullfile(dataDir, subID, 'raw', cntfn);
    fileID = [subID, '_', fnsuffix];
    if verbose
        disp("===================================================")
        disp([mHead, 'Current subject file ID: "', fileID, '"'])
        disp(" ")
    end
end

% Specify analysis output directory
outputDir = fullfile(dataDir, subID, 'analysis');
if ~exist(outputDir, 'dir')
    mkdir(outputDir) % make the directory if non-existant
end

end
