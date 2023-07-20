function [ edf_annot_channel, staging_offset_times ] = SleepEEG_buildannot( edfFN, staging_comment )
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - BUILDANNOT**
%
% - this is a function used to construct an annotation channel signal based
% on an EDF+ file exported from the Natus clinical SleepWorks system at
% the MGH Sleep Center. The purpose of this is to add a dummy comment every
% 5s so that they show up as annotations when the EDF is exported back
% into the Natus clinical SleepWorks system for us to extract the sleep
% stage scoring after the sleep tech has scored the imported study record
%
% This function automatically detects the number of samples in a single
% data record for the annotation channel and keeps it at that number. It
% will fill in the data records that don't have a comment in the original
% exported EDF+ file with dummy staging comments. As a result, The
% intervals between these dummy comments may not be exactly 5s, and will
% be shifted by the duration of a data record in order to maintain the same
% maximal samples in a single data record for the annotation channel.
%
% The output is a vector of signed integers of double type read as int16.
% You can treat this as the annotation channel data and feed it into
% blockEdfWrite() directly along with other channels that are also read as
% int16 (which is the default of SleepEEG_loadedf()).
%
% ===> Staging Comment abbreviated symbol inserted:
if nargin < 2 || isempty(staging_comment)
    staging_comment = 'Sleep Staging';
end
%
% Update on 02/04/2021: we are going to try with a new Staging Comment and
% see if we can let Natus recognize as a custom-made event type. This way
% when clinicians are scoring, they can filter for that particular event
% type and hide them, while preserving the ability to look at the other
% comments sleep techs noted over the course of a night.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - edfFN:        a string of filename of the exported EDF+ file.
%                           Usually should have _deidentified.edf suffix.
%
%           - staging_comment:
%                           a string to be used as inserted comments for
%                           exporting sleep stages after scoring.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - edf_annot_channel:
%                           a vector of annotation channel converted into
%                           the physical value scale as would be done by
%                           the SleepEEG_loadedf() function on a raw EDF+
%                           file data. This is kept at the same length as
%                           the signal in the SignalCell array as would be
%                           exported by calling SleepEEG_loadedf() with 4th
%                           outputs.
%
%           - staging_offset_times:
%                           a vector of double type containing the offset
%                           times of the Staging Comments inserted into the
%                           outputted edf_annot_channel. Offset is
%                           calculated as difference from the beginning of
%                           the recording in the EDF+ file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%%
%------------------------------------------------------------- Input check
% Check that first argument is a string
if ~ischar(edfFN)
    msg = ('Input argument is not a string.');
    error(msg);
end
%--------------------------------------------------- Load File Information
% Open file for reading
[fid, msg] = fopen(edfFN);
% Proceed if file is valid
if fid <0; error(msg); end

%%
%------------------------------------------------------------- Load Header
try
    % Load header information in one call
    edfHeaderSize = 256;
    [A, count] = fread(fid, edfHeaderSize);
    assert(count == edfHeaderSize, ['Not all ', num2str(edfHeaderSize) ' are successfully loaded!'])
catch exception
    msg = 'File load error. Cannot read header information. Check available memory.';
    error(msg);
end
%---------------------------------------------------- Process Header Block
% These specifications are hard-coded by the EDF format.
% See https://www.edfplus.info/specs/edf.html for descriptions.

% Create array/cells to create struct with loop
headerVariables = {...
    'edf_ver';            'patient_id';         'local_rec_id'; ...
    'recording_startdate';'recording_starttime';'num_header_bytes'; ...
    'reserve_1';          'num_data_records';   'data_record_duration';...
    'num_signals'};
headerVariablesConF = {...
    @strtrim;   @strtrim;   @strtrim; ...
    @strtrim;   @strtrim;   @str2num; ...
    @strtrim;   @str2num;   @str2num;...
    @str2num};
headerVariableSize = [8; 80; 80; 8; 8; 8; 44; 8; 8; 4];
headerVarLoc = vertcat(0, cumsum(headerVariableSize));
headerSize = sum(headerVariableSize);

% Create Header Structure
header = struct();
for h = 1:length(headerVariables)
    conF = headerVariablesConF{h};
    value = conF(char(A(headerVarLoc(h)+1:headerVarLoc(h+1))'));
    header.(headerVariables{h}) = value;
end

% End Header Load section

%%
%------------------------------------------------------ Load Signal Header
try
    % Load signal header into memory in one load
    edfSignalHeaderSize = header.num_header_bytes - headerSize;
    [A, count] = fread(fid, edfSignalHeaderSize);
    assert(count == edfSignalHeaderSize, ['Not all ', num2str(edfSignalHeaderSize) ' are successfully loaded!'])
catch exception
    msg = 'File load error. Cannot read signal header information. Check available memory.';
    error(msg);
end
%--------------------------------------------- Process Signal Header Block
% These specifications are hard-coded by the EDF format.
% See https://www.edfplus.info/specs/edf.html for descriptions.

% Create array/cells to create struct with loop
signalHeaderVar = {...
    'signal_labels'; 'tranducer_type'; 'physical_dimension'; ...
    'physical_min'; 'physical_max'; 'digital_min'; ...
    'digital_max'; 'prefiltering'; 'samples_in_record'; ...
    'reserve_2' };
signalHeaderVarConvF = {...
    @strtrim; @strtrim; @strtrim; ...
    @str2num; @str2num; @str2num; ...
    @str2num; @strtrim; @str2num; ...
    @strtrim };
num_signal_header_vars = length(signalHeaderVar);
num_signals = header.num_signals;
signalHeaderVarSize = [16; 80; 8; 8; 8; 8; 8; 80; 8; 32];
signalHeaderBlockSize = sum(signalHeaderVarSize)*num_signals; %#ok<NASGU>
signalHeaderVarLoc = vertcat(0, cumsum(signalHeaderVarSize*num_signals));
signalHeaderRecordSize = sum(signalHeaderVarSize); %#ok<NASGU>

% Create Signal Header Struct
signalHeader = struct(...
    'signal_labels', {},'tranducer_type', {},'physical_dimension', {}, ...
    'physical_min', {},'physical_max', {},'digital_min', {},...
    'digital_max', {},'prefiltering', {},'samples_in_record', {},...
    'reserve_2', {});

% Get each signal header variable
for v = 1:num_signal_header_vars
    varBlock = A(signalHeaderVarLoc(v)+1:signalHeaderVarLoc(v+1))';
    varSize = signalHeaderVarSize(v);
    conF = signalHeaderVarConvF{v};
    for s = 1:num_signals
        varStart = varSize*(s-1)+1;
        varEnd = varSize*s;
        value = conF(char(varBlock(varStart:varEnd)));
        signalHeader(s).(signalHeaderVar{v}) = value;
    end
end

% End Signal Header Load Section

%%
%--------------------------------------------- Prepare to Load Annotations
% Setting up size and indexing variables
annot = arrayfun(...
    @(x)strcmp(signalHeader(x).signal_labels,'EDF Annotations'),...
    1:header.num_signals, 'UniformOutput', false);
annotidx = find(cell2mat(annot)==1);

edfSignalSizes = arrayfun(...
    @(x)signalHeader(x).samples_in_record, 1:header.num_signals);
edfRecordSize = sum(edfSignalSizes);
%--------------------------------------------- Process Annotations channel
% Rewind the fid pointer to the start of the file
frewind(fid);
% Burn out the header lines
[~, ~] = fread(fid, header.num_header_bytes);
% find the offset time of the Lights Off comment and the last comment
lightoff_entry = [];
lastcomment_entry = [];
for r = 1:header.num_data_records
    % first burn everything before the annotation channel
    [~, ~] = fread(fid, sum(edfSignalSizes(1:annotidx-1)), 'int16');
    % read in the annotation channel
    [a, count] = fread(fid, edfSignalSizes(annotidx)*2);
    assert(count == edfSignalSizes(annotidx)*2, ['Not all ', num2str(edfSignalSizes(annotidx)*2) ' are successfully loaded!'])
    % this should have exhausted one data record
    assert(edfRecordSize - sum(edfSignalSizes(1:annotidx)) == 0, 'Non-empty remaining samples in a data record.')
    % find the data record containing the Lights Off comment
    tempstr = char(a');
    if contains(tempstr, 'Lights Off') % Lights Off comment
        lightoff_entry = tempstr;
    elseif r == header.num_data_records % last comment
        lastcomment_entry = tempstr;
    end
end
assert(~isempty(lightoff_entry), 'No Lights Off comment was found in the exported EDF+ file. Please check in the Natus system!')
assert(contains(lastcomment_entry, char(43)), 'Last comment does not have an offset time. Please check in the annotation channel!')

% Now parse the annotation strings into Time-stamped Annotations Lists (TALs)
% and find the Lights Off offset time. We also need to find the last offset
% time for comments.
lightoff_time = [];
tempstr = strsplit(lightoff_entry, char([20, 0]));
for ii = 1:length(tempstr)
    currentstr = tempstr{ii};
    if contains(currentstr, 'Lights Off')
        parselist = strsplit(currentstr, {char(43), char(20)});
        onsetstr = parselist{find(~cellfun(@isempty, parselist), 1, 'first')};
        assert(~contains(onsetstr, char(21)), 'The Lights Off comment has a non-zero duration. Something is wrong!')
        lightoff_time = str2num(onsetstr);
        %also need to get the number of decimals in the Lights Off offset
        %time
        lightoff_split = strsplit(onsetstr, '.');
        % Update 06/18/2021: if Lights Off was manually moved after placing
        % by sleep tech, it will have a small decimal place. We force it to
        % be 6 at minimal.
        if length(lightoff_split) == 1
            lightoff_split{2} = '000000';
        end
        numofdeci = max(6,length(lightoff_split{2}));
    end
end
assert(~isempty(lightoff_time), 'Lights Off offset time was not successfully extracted.')

lastcomment_time = [];
tempstr = strsplit(lastcomment_entry, char([20, 0]));
for ii = 1:length(tempstr)
    currentstr = tempstr{ii};
    if contains(currentstr, char(43))
        parselist = strsplit(currentstr, {char(43), char(20)});
        onsetstr = parselist{find(~cellfun(@isempty, parselist), 1, 'first')};
        assert(~contains(onsetstr, char(21)), 'The last comment has a non-zero duration. Something is wrong!')
        lastcomment_time = str2num(onsetstr);
    end
end
assert(~isempty(lastcomment_time), 'Last comment offset time was not successfully extracted.')

% Now we construct a series of offset timepoints for our dummy staging
% comments.
% based on the lightoff_time, we need to figure out when should Staging
% Comments begin. We will roll back to the beginning of the start of the
% 30-s epoch that contains this Lights Off comment
assert(numofdeci > 5, 'Too few decimal places for the lightoff_time, please check!')
start_comment_time = double(idivide(lightoff_time, int16(30))) * 30 - 30 +...
    str2num(lightoff_split{2}(end-1:end))/10^numofdeci; % add a small decimal shift
offset_seq = start_comment_time:5:lastcomment_time; % insert every 5 s
assert(offset_seq(1) == start_comment_time && ...
    offset_seq(end) > lastcomment_time-5, 'Boundaries of offset timepoints are incorrect.')
noffset = length(offset_seq);

%%
%------------------------------------------- Construct Annotations channel
% Rewind the fid pointer to the start of the file
frewind(fid);
% Burn out the header lines
[~, ~] = fread(fid, header.num_header_bytes);
% Inserting dummy staging comments into the right places
n_inserted = 0;
annotation_channel = zeros(signalHeader(annotidx).samples_in_record * header.num_data_records*2, 1);
for r = 1:header.num_data_records
    % first burn everything before the annotation channel
    [~, ~] = fread(fid, sum(edfSignalSizes(1:annotidx-1)), 'int16');
    % read in the annotation channel
    [a, count] = fread(fid, edfSignalSizes(annotidx)*2);
    assert(count == edfSignalSizes(annotidx)*2, ['Not all ', num2str(edfSignalSizes(annotidx)*2) ' are successfully loaded!'])
    % this should have exhausted one data record
    assert(edfRecordSize - sum(edfSignalSizes(1:annotidx)) == 0, 'Non-empty remaining samples in a data record.')
    % convert to character from double strings
    current_char = char(a');
    % remove annotations to be ignored when importing
    current_char = ignore_TALs(current_char, {'Oxygen Desaturation', 'Observation Note'});
    new_char = current_char;
    % find the offset time corresponding to the current data record
    split_list = strsplit(current_char, char([20, 0]));
    assert(contains(split_list{1}, {char(43), char(20)}), 'First position TAL does not correspond to an offset time')
    tempstr = strsplit(split_list{1}, {char(43), char(20)});
    current_offset = str2num(tempstr{find(~cellfun(@isempty, tempstr), 1, 'first')});
    assert(current_offset > 0, 'Not a valid offset time value!')
    % now compare current_offset with the next offset_seq on the list to
    % decide whether to inset a Staging Comment at this data record
    if n_inserted < noffset && current_offset > offset_seq(n_inserted+1)-0.5
        % we will insert it here, but the issue is there might be a comment
        % here already. We have to find the right place to insert.
        current_last_comment = split_list{end-1};
        % find the position of the last comment in this data record
        last_start_index = strfind(current_char, current_last_comment);
        % extract the offset time of the last comment
        last_comment_split = strsplit(current_last_comment, {char(43), char(20)});
        last_comment_offset_str = last_comment_split{find(~cellfun(@isempty, last_comment_split), 1, 'first')};
        if contains(last_comment_offset_str, char(21))
            tempsplit = strsplit(last_comment_offset_str, char(21));
            last_comment_offset_str = tempsplit{1};
        end
        last_comment_offset = str2num(last_comment_offset_str);
        % modify the offset value in offset_seq if last comment is larger.
        % We manually make sure the offset value of our Staging Comment has
        % the largest offset within a data record
        if last_comment_offset > offset_seq(n_inserted+1)
            offset_seq(n_inserted+1) = last_comment_offset + (current_offset+0.5 - last_comment_offset)/2;
        end
        % we start replacing two more characters away because we have
        % removed the char([20, 0]) during strsplit.
        insert_comment = [char(43), eval([char("sprintf('%."), char(sprintf("%df', offset_seq(n_inserted+1))", numofdeci))]),...
            char(21), '1', char(20), staging_comment, char([20, 0])]; %  Update on 02/04/2021: give it a 1s duration
        start_index = last_start_index+length(current_last_comment)+2;
        end_index = start_index+length(insert_comment)-1;
        new_char(start_index:end_index) = insert_comment;
        % check that the maximum length is not exceeded
        assert(length(new_char) == length(current_char), 'Maximum sample size exceeded during inserting comments!')
        % update the number of inserted comments
        n_inserted = n_inserted + 1;
    end
    new_a = double(new_char');
    annotation_channel((r-1)*edfSignalSizes(annotidx)*2+1:r*edfSignalSizes(annotidx)*2) = new_a;
end

% Close file explicitly
if fid > 0; fclose(fid); end

%%
% Quality checks of the constructed Annotations channel
% 1) basic check of signal length
assert(n_inserted == noffset, 'Number of inserted comments does not match what should have been inserted.')
assert(length(annotation_channel)/2 == signalHeader(end).samples_in_record * header.num_data_records,...
    'Length of the annotation channel does not match that specified by header and signalHeader.')
if isfile('Decoding_test_delete_later.edf'); delete 'Decoding_test_delete_later.edf'; end
[fid, msg] = fopen('Decoding_test_delete_later.edf', 'w+');
if fid <0; error(msg); end
fwrite(fid, int8(annotation_channel), 'int8');
fclose(fid);
[fid, msg] = fopen('Decoding_test_delete_later.edf', 'r+');
if fid <0; error(msg); end
[edf_annot_channel, ~] = fread(fid, 'int16');
assert(length(edf_annot_channel) == signalHeader(end).samples_in_record * header.num_data_records,...
    'Length of the annotation channel when written into EDF+ file does not match that specified by header and signalHeader.')

% 2) decode the annotation channel
frewind(fid);
annotation = {};
for r = 1:header.num_data_records
    [a, count] = fread(fid, edfSignalSizes(annotidx)*2);
    assert(count == edfSignalSizes(annotidx)*2, ['Not all ', num2str(edfSignalSizes(annotidx)*2) ' are successfully loaded!'])
    % this should have exhausted one data record
    assert(edfRecordSize - sum(edfSignalSizes(1:annotidx)) == 0, 'Non-empty remaining samples in a data record.')
    % now let's parse the annotation strings into Time-stamped
    % Annotations Lists (TALs)
    tempstr = char(a');
    if contains(tempstr, char([20, 0]))
        tempstr = strsplit(tempstr, char([20, 0]));
        annotation(size(annotation,1)+1:size(annotation,1)+length(tempstr), 1) = tempstr';
    else
        error('Current annotation entry in the data record does not have char([20, 0])')
    end
end
annotation = SleepEEG_annot(header, annotation);

% 3) check the ordering of all comment offset times
all_offset_times = zeros(1,length(annotation));
deleteidx = [];
for ii = 1:length(annotation)
    all_offset_times(ii) = annotation(ii).offset;
    if ~strcmp(annotation(ii).annotation, staging_comment)
        deleteidx = [deleteidx, ii]; %#ok<*AGROW>
    end
end
assert(all(all_offset_times == sort(all_offset_times)), 'Ordering of offset times is incorrect.')

% 4) check the intervals of inserted Staging Comments
annotation(deleteidx) = [];
assert(length(annotation) == noffset, 'Number of Staging Comments inserted is incorrect.')
staging_offset_times = zeros(1,length(annotation));
for ii = 1:length(annotation)
    staging_offset_times(ii) = annotation(ii).offset;
end
assert(all(diff(staging_offset_times) > 5-header.data_record_duration & diff(staging_offset_times) < 5+header.data_record_duration),...
    'Intervals between Staging Comments seem to be incorrect.')

% all tests passed, delete temporary file
delete 'Decoding_test_delete_later.edf'

%% Convert from digital to physical values like SleepEEG_loadedf()
% we add this section so that blockEdfWrite can properly handle the
% conversion back to digital when writing all signals into an EDF+ file.

% Get scaling factors
dig_min = signalHeader(annotidx).digital_min;
dig_max = signalHeader(annotidx).digital_max;
phy_min = signalHeader(annotidx).physical_min;
phy_max = signalHeader(annotidx).physical_max;

% Convert from digital to physical values
value = (edf_annot_channel-dig_min)/(dig_max-dig_min);
value = value.*double(phy_max-phy_min)+phy_min;

% Update the output vector
edf_annot_channel = value;

end

%% HELPER FUNCTIONS
function [ current_char ] = ignore_TALs(current_char, ignored_annotation)
% added on 02/04/2021: we add a functionality here to remove all
% the 'Oxygen Desaturation' and 'Observation Note' comments to help
% with running Batch Analyzer on the imported EDF file
if nargin < 2
    ignored_annotation = {'Oxygen Desaturation', 'Observation Note'};
end

for ignore_string = ignored_annotation
    ignore_char = ignore_string{1};
    if contains(current_char, ignore_char)
        % split the character into TALs and check offset time validity
        split_list = strsplit(current_char, char([20, 0]));
        assert(contains(split_list{1}, {char(43), char(20)}), 'First position TAL does not correspond to an offset time')
        tempstr = strsplit(split_list{1}, {char(43), char(20)});
        current_offset = str2num(tempstr{find(~cellfun(@isempty, tempstr), 1, 'first')});
        assert(current_offset > 0, 'Not a valid offset time value!')
        
        % find out which TALs contains the comment to remove
        remove_index = 0;
        for ii = 1:length(split_list)
            currentstr = split_list{ii};
            if contains(currentstr, ignore_char)
                remove_index = ii;
                comment_to_remove = currentstr; % get the comment
            end
        end
        assert(remove_index ~= 0, 'Could not find the TAL to remove.')
        % find out the index at which the TAL starts
        remove_start_index = strfind(current_char, comment_to_remove);
        % we need to replace two more indices because at the end of each
        % TAL there are the char([20, 0]) that got removed during strsplit
        replace_comment = char(zeros(1, length(comment_to_remove)+2));
        % update current_char
        current_char(remove_start_index:remove_start_index+length(replace_comment)-1) = [];
        current_char = [current_char, replace_comment];
    end
end
end
