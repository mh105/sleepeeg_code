function [ channelNum, srate ] = SleepEEG_extractChan(subID, fnsuffix, channelNum, project, dataDir, datafn, fileID)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - EXTRACTCHAN**
%
% - function used to extract several channels from overnight sleep EEG
% data recorded with ANT systems (asalab or eego lab)
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
%           - channelNum:   a vector containing channel numbers to extract,
%                           e.g. [25, 84].
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
% ##### Outputs:
%
%           - channelNum:   a vector containing channel numbers to extract,
%                           repeating so that reconChan can take as input.
%
%           - srate:        sampling rate of the extracted EEG channel, of
%                           double type
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 3
    channelNum = [25, 84]; % by default extract Z10-25, Z2-84
    project = '';
elseif nargin < 4
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('--------------------------------')
disp([mHead 'SleepEEG_extractChan()']);
disp('--------------------------------')

%% Report extracted channels
disp([mHead, 'Extracting channels in original sampling frequency: ', num2str(channelNum)])

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
disp([mHead, 'Loading and Extracting Channels...'])
disp(' ')

% Create a temporary directory to hold the segments data
tempDir = fullfile(dataDir, subID, 'temp');
if ~exist(tempDir, 'dir')
    mkdir(tempDir) % make the directory if non-existant
end
disp([mHead, 'Temporary folder to dump channel pieces:'])
fprintf(mSpace)
disp(tempDir)

% construct 30-min segment sample points for looping
sample_point = 0:cnt_info.sample_rate*60*30:cnt_info.sample_count;
sample_point = [sample_point cnt_info.sample_count]; % add in the last sample point
disp(' ')
if length(sample_point) > 2
    disp([mHead, 'A total of ' num2str(length(sample_point)-1) ' segments of 30min'])
else
    disp([mHead, 'A total of ' num2str(length(sample_point)-1) ' segments of ' num2str(record_time) 'min'])
end

tic
for i = 1:length(sample_point)-1
    % load in one segment of data
    if i == 1
        sample1 = 1;
    else
        sample1 = sample_point(i)+1;
    end
    sample2 = sample_point(i+1);
    disp([mHead, 'Loading the ' num2str(i) 'th segment...'])

    % use pop_loadeep_v4.m function to load in the segment
    EEG = pop_loadeep_v4(datafn, 'sample1', sample1, 'sample2', sample2);

    % store the sampling rate
    srate = EEG.srate;

    for j = channelNum % loop through each channel number
        chan_data = EEG.data(j, :);

        % save the segment for ONE channel
        filename = [fileID '_Channel_' num2str(j) '_' num2str(i) '_segment'];
        save(fullfile(dataDir, subID, 'temp', filename), 'chan_data')
    end

    % clearvars to save memory
    clearvars EEG chan_data
end

disp([mHead, 'Total time taken in Loading and Extracting channels...'])
fprintf(mSpace)
toc

end
