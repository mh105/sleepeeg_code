function [ channelfn ] = SleepEEG_reconChan(subID, fnsuffix, channelNum, project, dataDir, datafn, fileID)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - RECONCHAN**
%
% - function used to extract one occipital channel from overnight sleep EEG
% data recorded with ANT systems
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep".
%
%                           ***the final file name is in the form: subID_fnsuffix.cnt***
%
%           - channelNum:   a vector containing channel numbers to extract,
%                           e.g. [25, 84]. Required.
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
%           - dataDir:      directory path containing all subjects' data.
%
%           - datafn:       full path to the .cnt file including preceding
%                           directories. Can be obtained from
%                           SleepEEG_configDir(subID, fnsuffix, verbose)
%
%           - fileID:       name of the .cnt file (no .cnt suffix).
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - channelfn:    a cell array of filenames of saved channel
%                           data.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

assert(exist('channelNum', 'var')==1, 'No channel number provided for reconstruction!')

if nargin < 4
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('------------------------------')
disp([mHead 'SleepEEG_reconChan()']);
disp('------------------------------')

%% Report reconstructed channels
disp([mHead, 'Reconstructing channels in original sampling frequency: ', num2str(channelNum)])

%% Define Directories of Codes and Data Folders
if ~exist('dataDir', 'var')
    [ dataDir, datafn, fileID, ~ ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

if ~exist('fileID', 'var')
    [ ~, ~, fileID ] = SleepEEG_configDir(subID, fnsuffix, false, project);
end

%% Report EEG Record Length
cnt_info = eepv4_read_info(datafn);
tot_sample_point = cnt_info.sample_count;
record_time = cnt_info.sample_count/cnt_info.sample_rate/60;

disp(' ')
disp([mHead, 'Total recording duration = ' num2str(record_time) ' minutes, ', num2str(record_time/60) ' hours.'])
disp([mHead, 'Total sample points = ' num2str(tot_sample_point) '.'])

%% Load Data In Pieces
disp(' ')
disp([mHead, 'Loading and Concatenating Channels...'])
disp(' ')

% Create a channel directory to store single channel data
chanDir = fullfile(dataDir, subID, 'channel');
if ~exist(chanDir, 'dir')
    mkdir(chanDir) % make the directory if non-existant
end
disp([mHead, 'Folder to store single channel data:'])
fprintf(mSpace)
disp(chanDir)

% construct 30-min segment sample points for looping
sample_point = 0:cnt_info.sample_rate*60*30:cnt_info.sample_count;
sample_point = [sample_point cnt_info.sample_count]; % add in the last sample point
disp(' ')
if length(sample_point) > 2
    disp([mHead, 'A total of ' num2str(length(sample_point)-1) ' segments of 30min'])
else
    disp([mHead, 'A total of ' num2str(length(sample_point)-1) ' segments of ' num2str(record_time) 'min'])
end

% Create a cell for storing the saved channel data file names
channelfn = cell(1, length(channelNum));

tic
for j = channelNum % loop through each channel number
    % Initialize the channel data vector
    channel_data = zeros(1, cnt_info.sample_count);

    % load each segment sequentially and save as one file
    for i = 1:length(sample_point)-1

        % load in one segment of data
        filename = [fileID '_Channel_' num2str(j) '_' num2str(i) '_segment'];
        load(fullfile(dataDir, subID, 'temp', filename), 'chan_data')
        disp([mHead, 'Loading the ' num2str(i) 'th channel segment...'])

        % specify indices in channel data
        if i == 1
            sample1 = 1;
        else
            sample1 = sample_point(i)+1;
        end
        sample2 = sample_point(i+1);

        channel_data(sample1:sample2) = chan_data;
        clearvars chan_data

    end

    % save the segment for ONE channel
    filename = [fileID '_Channel_' num2str(j) '_all_segment'];
    save(fullfile(dataDir, subID, 'channel', filename), 'channel_data')

    % add file name to the cell array for returning function output
    channelfn{channelNum==j} = fullfile(dataDir, subID, 'channel', filename);
end

% Clear the temp folder to save disk space
which_dir = fullfile(dataDir, subID, 'temp');
dinfo = dir(which_dir);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(which_dir, {dinfo.name});
delete(filenames{:})

disp([mHead, 'Total time taken in Reconstructing channels...'])
fprintf(mSpace)
toc

end
