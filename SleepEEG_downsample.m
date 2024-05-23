function [ EEG ] = SleepEEG_downsample(subID, fnsuffix, dsrate, verbose, oversave, dataDir, datafn)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - DOWNSAMPLE**
%
% - loads a recording from a .cnt file and downsamples the eego mylab high
% density EEG recording to 'dsrate' Hz in preparation for alignment with
% clinical PSG recordings
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
%           - dsrate:       whether to downsample and desired downsampling rate.
%                           default: [true, 500]
%
%           - verbose:      whether print messages during processing.
%                           default: true
%
%           - oversave:     whether to overwrite existing .set files and
%                           save the EEG structure again onto disk
%                           default: false
%
%           - dataDir:      directory path containing all subjects' data.
%
%           - datafn:       full path to the .cnt file including preceding
%                           directories. Can be obtained from
%                           SleepEEG_configDir(subID, fnsuffix, verbose)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - EEG:          an EEGLAB structure containing all information
%                           of the recording in .cnt file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    dsrate = [true, 500];
    verbose = true;
    oversave = false;
elseif nargin < 4
    verbose = true;
    oversave = false;
end

if ~exist('dataDir', 'var')
    [dataDir, datafn, fileID, ~] = SleepEEG_configDir(subID, fnsuffix, verbose);
end

if ~exist('fileID', 'var')
    fileID = [subID, '_', fnsuffix];
end

% If the path to ANT_interface_code was not added, add it now
if ~exist('ANT_interface_readcnt', 'file') == 2
    SleepEEG_addpath(matlabroot);
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('-------------------------------')
disp([mHead 'SleepEEG_downsample()']);
disp('-------------------------------')

%% Configure filename and filepath
filename = [fileID, '.cnt'];
pathsplit = strsplit(datafn, filename);
filepath = pathsplit{1};

subfolder = 'set';
savefilepath = fullfile(dataDir, subID, subfolder);

if ~exist(savefilepath, 'dir')
    mkdir(savefilepath) % make the directory if non-existant
end
if verbose
    disp('Destination folder to store downsampled data:')
    disp(' ')
    disp(savefilepath)
    disp(' ')
end

if dsrate(1) % if downsampling flag is true
    savefn = [subID, '_', fnsuffix, '_ds', num2str(dsrate(2)), '_Z3']; % Specify .set file name
else % save at original sampling rate
    savefn = [subID, '_', fnsuffix, '_Z3']; % Specify .set file name
end

%% Read .cnt file with downsampling

% we skip reading .cnt and downsampling if the desired file is already on
% disk and oversave is set to be false
if ~oversave && isfile(fullfile(dataDir, subID, subfolder, [savefn, '.set']))
    if verbose
        disp([mHead, savefn, '.set already exists, and oversave is false.']);
        disp([mHead, 'Reading directly from disk, not overwriting.']);
        disp(' ')
    end
    
    EEG = ANT_interface_loadset([savefn, '.set'], fullfile(dataDir, subID, subfolder), verbose, false);
    
elseif oversave && isfile(fullfile(dataDir, subID, subfolder, [savefn, '.set']))
    if verbose
        disp([mHead, savefn, '.set already exists, but oversave is true.']);
        disp([mHead, 'Reading raw .cnt file now.']);
        disp(' ')
    end
    
    EEG = ANT_interface_readcnt(filename, filepath, dsrate, verbose);
    
    % Save downsampled data
    if verbose; disp([mHead, 'Oversave is true, overwriting ', savefn, '.set']); end
    EEG = ANT_interface_saveset(EEG, savefn, savefilepath, verbose);
    
else % otherwise we need to use ANT_interface_readcnt() to load the data
    if verbose
        disp([mHead, savefn, '.set is not found, reading raw .cnt file now.']);
        disp(' ')
    end
    
    EEG = ANT_interface_readcnt(filename, filepath, dsrate, verbose);

    % Save downsampled data
    if verbose; disp([mHead, savefn, '.set does not exist, saving data now.']); end
    EEG = ANT_interface_saveset(EEG, savefn, savefilepath, verbose);
    
end

end
