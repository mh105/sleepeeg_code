function [] = SleepEEG_convertscoring(txtFN, edfFN, alignedfFN, filepath, rvsFN, staging_comment)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - CONVERTSCORING**
%
% - this is a function used to read in the exported scoring text file after
% importing an EDF/EDF+ study back to the Natus system. After scoring is
% done by the sleep techs, the entire comment section is exported as a .txt
% file, and each comment is labelled with an epoch number and a scored
% sleep stage. We will read in the text file and update the Time column to
% be off-sets from the beginning of HD-EEG recordings. This converted
% scoring text file will then be imported with the HD-EEG data in BIDS
% format and stored as a MNE raw object. Various quality checks are
% completed to ensure accuracy of the extracted scoring
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - txtFN:        filename of the scoring .txt file.
%
%           - edfFN:        filename of the exported .edf file that ends
%                           with _deidentified.edf.
%
%           - alignedfFN:   filename of the aligned .edf file that ends
%                           with _aligned.edf.
%
%           - filepath:     path to the "clinical" folder under the
%                           subject's folder.
%
%           - rvsFN:        filename of the reverse alignment .mat file
%                           ouputted from trigger alignment code.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 6 || isempty(staging_comment)
    staging_comment = 'Sleep Staging';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead)); %#ok<NASGU>

%%
disp('-----------------------------------')
disp([mHead 'SleepEEG_convertscoring()']);
disp('-----------------------------------')

%% Read the text file
disp([mHead, 'Loading the raw comment file...'])
% temporarily disable the UTF-16LE warning message
warning('off', 'MATLAB:iofun:UnsupportedEncoding')

% Load the comments file (.txt) into a cell array
delimiter = '\t';
startRow = 7;
formatSpec = '%q%q%q%q%q%[^\n\r]';
fileID = fopen(fullfile(filepath, txtFN),'r','n','UTF-16LE');
fseek(fileID, 2, 'bof');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
commentfile = table(dataArray{1:end-1}, 'VariableNames', {'Time','Epoch','Stage','Duration','Title'});
clearvars delimiter startRow formatSpec fileID dataArray ans;

% turn UTF-16LE warning message back on
warning('on', 'MATLAB:iofun:UnsupportedEncoding')

% Find the staging comments and recording start time
stagingidx = false(size(commentfile,1), 1);
for ii = 1:size(commentfile,1)
    if strcmp(commentfile.Title{ii}, 'Beginning of Recording')
        record_starttime = real2passed(commentfile.Time{ii});
    end
    if strcmp(commentfile.Title{ii}, staging_comment)
        stagingidx(ii) = true;
    end
end

%% Process all comments
disp([mHead, 'Extracting Scoring Comments...'])
offset_list = arrayfun(@(x) real2passed(x{1}, record_starttime), commentfile.Time);
epoch_list = arrayfun(@(x) getepoch(x{1}), commentfile.Epoch);
stage_list = arrayfun(@(x) getstage(x{1}), commentfile.Stage);

% sanity check that staging comment duration should be 0:01.0
if strcmp(staging_comment, 'Sleep Staging')
    correctduration = arrayfun(@(x) strcmp(x{1}, '0:01.0'), commentfile.Duration(stagingidx));
    assert(all(correctduration(1:end-1)==1), 'Unexpected duration detected for Staging Comments.')
end

%% Quality checks of extracted staging
disp([mHead, 'Checking the quality of extracted staging scoring...'])
% 1. Check that the Staging Comments are spaced roughly 5s apart
interval_list = diff(offset_list(stagingidx));
assert(all(interval_list>=4 & interval_list<=6), 'Staging Comments are not 5s apart.')
% 2. Check that we have exactly 6 Staging Comments for every epoch
temp_epoch_list = epoch_list(stagingidx); % except for the 1st and last epochs
temp_epoch_list(epoch_list(stagingidx) == 0 | epoch_list(stagingidx) == max(epoch_list(stagingidx))) = [];
uv = unique(temp_epoch_list);
epoch_count = histcounts(temp_epoch_list, [uv; uv(end)+1]);
assert(all(epoch_count(uv~=0)==6), 'Not exactly 6 Staging Comments for every epoch after Lights Off.')
% 3. Check that we have the same number of Staging Comments as inserted
[ ~, staging_offset_times ] = SleepEEG_buildannot(fullfile(filepath, edfFN));
assert(length(staging_offset_times) == length(stage_list(stagingidx)), 'Different numbers of Staging Comments exported from inserted.')
% 4. Check that the offset times extracted match the more precisely inserted
offset_diff = offset_list(stagingidx) - staging_offset_times'; % we allow < 0.5s discrepancy
assert(all(offset_diff < 0.5), 'Extracted offset times are inaccurate compared to inserted times.')

%% Reverse align offset times to the HD-EEG timeframe
disp([mHead, 'Reverse alignment to the HD-EEG timeframe...'])
% Load the reverse alignment structure saved in trigger alignment code
load(fullfile(filepath, rvsFN), 'rvsalign_store')

% In some of the early recordings, stochastic triggers were not
% available, so Amanda had to perform manual alignment that
% unfortunately resulted in incorrect indices saved in rvsalign_store.
% In these cases, we need to search for the matching points using the
% C3 channel.

% Load the C3 channel in the aligned EDF+ file
[header, signalHeader, signalCell] = SleepEEG_loadedf(fullfile(filepath, alignedfFN), {'C3'});
Fs = signalHeader.samples_in_record / header.data_record_duration;
assert(Fs == 500, 'Sampling frequency is not 500 Hz.')
aligned_signal = signalCell{1};

if rvsalign_store.truncate_endidx > rvsalign_store.original_length
    % truncate_endidx > original_length is the indication that this
    % recording was manually aligned because in the code used for
    % manual alignment, original_length = truncate_endidx - a number.
    % This could never happen in automatic alignment because
    % truncate_endidx is supposed to cut the original EEG data
    disp([mHead, '[!] Manual alignment detected, need to search for matching points...'])
    
    % remove the padding on both ends of aligned signal in edf file
    aligned_signal_nopadding = aligned_signal(rvsalign_store.begin_pad+1:length(aligned_signal)-rvsalign_store.end_pad);
    
    % now we need to load the set file and grab the EEG channel used to
    % emulate C3 in the aligned edf - this takes a while to load.
    subID = txtFN(1:strfind(txtFN,'_night')-1);
    night = txtFN(strfind(txtFN,'_night')+6);
    filename = [subID, '_night' night '_Sleep_ds500_Z3.set'];
    set_filepath = strrep(filepath, 'clinical', 'set');
    EEG_new = ANT_interface_loadset(filename, set_filepath, true, false);
    originalEEG = EEG_new.data(label2num('LA2',EEG_new.chanlocs),:); % C3 corresponds to LA2
    
    % find where the aligned_signal_nopadding starts and ends
    search_tic = tic;
    length_diff = length(originalEEG) - length(aligned_signal_nopadding);
    assert(length_diff > 0, 'truncated EEG is longer than original EEG. Something is wrong.')
    all_diff_sum = [];
    for ii = 1:length_diff
        all_diff_sum(ii) = sum(abs(originalEEG(ii:ii+length(aligned_signal_nopadding)-1)-aligned_signal_nopadding')); %#ok<*AGROW>
    end
    [~, aligned_startidx] = min(all_diff_sum);
    disp([mHead, 'Time taken in searching for matching points:'])
    toc(search_tic)
    
    % figure out the offset difference between timeframes in seconds
    offset_adjustment = (aligned_startidx - (rvsalign_store.begin_pad+1)) / Fs;
    
else % with correctly behaving stochastic triggers, the rvsalign_store carries everything we need
    % figure out the cutting indices
    cut_begin = rvsalign_store.begin_pad - rvsalign_store.truncate_startidx + 2;
    cut_end = length(aligned_signal) - (rvsalign_store.end_pad - (rvsalign_store.original_length - rvsalign_store.truncate_endidx));
    
    % sanity check on the length of EEGvec_stage_channel
    assert(length(cut_begin:cut_end) == rvsalign_store.original_length, 'Length of EEG channel computed from cut_begin/end is incorrect.')
    
    % figure out the offset difference between timeframes in seconds
    offset_adjustment = -cut_begin / Fs;
    
end

%% Write an updated scoring text file
raw_scoring_fID = fopen(fullfile(filepath, txtFN),'r','n','UTF-16LE');
set_scoring_fn = fullfile(filepath, strrep(txtFN, '_scoring', '_scoring_set'));
set_scoring_fID = fopen(set_scoring_fn, 'wt');

% write the first 6 header lines from the original scoring text file
for ii = 1:6
    tline = fgetl(raw_scoring_fID);
    fprintf(set_scoring_fID, '%s\n', tline);
end
fclose(raw_scoring_fID);
fclose(set_scoring_fID);

% write the table with updated offset times in the Time column
commentfile.Time = offset_list + offset_adjustment;
writetable(commentfile, set_scoring_fn, 'Delimiter','\t', 'WriteVariableNames',false, 'Writemode','Append')

end

%% HELPER FUNCTIONS
function [ passedtime ] = real2passed(Time, record_starttime)
timesplit = strsplit(Time, ':');
if str2num(timesplit{1}) < 12 %#ok<*ST2NM>
    month = 1;day = 2;year = 1111; % pseudo M/D/Y for computing offsets
else
    month = 1;day = 1;year = 1111; % pseudo M/D/Y for computing offsets
end
current_time = datetime(year, month, day, str2num(timesplit{1}), str2num(timesplit{2}), str2num(timesplit{3}));
if nargin < 2
    passedtime = current_time;
else
    passedtime = seconds(current_time - record_starttime);
end
end

function [ epoch_num ] = getepoch(epochstr)
epoch_num = str2num(epochstr);
if isempty(epoch_num)
    epoch_num = 0;
end
end

function [ stagenum ] = getstage(stagestr)
% fill in unstaged string for early time points
if isempty(stagestr)
    stagestr = 'U';
end

% re-code the sleep stage to numbers per hypnoplot() function
% stages: 1 x S vector of stage values (0:Undefined, 5: Wake, 4:REM, 3:N1, 2:N2, 1:N3)
if strcmp(stagestr, 'W')
    stagenum = 5;
elseif strcmp(stagestr, 'R')
    stagenum = 4;
elseif strcmp(stagestr, 'N1')
    stagenum = 3;
elseif strcmp(stagestr, 'N2')
    stagenum = 2;
elseif strcmp(stagestr, 'N3')
    stagenum = 1;
elseif strcmp(stagestr, 'U')
    stagenum = 0;
else
    error('Sleep stage code includes non-decodable strings. Please check!')
end

end
