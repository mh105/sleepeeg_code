function [ new_annotation ] = SleepEEG_annot(header, annotation)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - ANNOT**
%
% - this is a function used to organize the annotation structure extracted
% from Natus clinical sleep EEG system using the SleepEEG_loadedf()
% function
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - header:       the header extracted by SleepEEG_loadedf().
%
%           - annotation:   the annotation structure extracted by
%                           SleepEEG_loadedf() but entries are in string
%                           format so inconvenient to use.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - new_annotation:
%                           an annotation structure containing the offset
%                           duration, and annotation comments, with fields
%                           cleaned up and organized.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%%
% Parse date/time of the start of recording
startdate = header.recording_startdate;
startdate = strsplit(startdate, '.');
if length(startdate) > 1
    month = str2num(startdate{1});
    day = str2num(startdate{2});
    year = str2num(startdate{3}) + 2000; % assuming data collected after year 2000
else
    month = 1;
    day = 1;
    year = 1111;
end
starttime = header.recording_starttime;
starttime = strsplit(starttime, '.');
hour = str2num(starttime{1});
minute = str2num(starttime{2});
second = str2num(starttime{3});

onset_t = datetime(year, month, day, hour, minute, second);

%%
% Clean up the annotation structure for uncessary cells
deleteidx = [];
for i = 1:length(annotation)
    currentstr = annotation{i};
    if ~contains(currentstr, char(43))
        deleteidx = [deleteidx, i];
    end
end
annotation(deleteidx) = [];

%%
% Now we will organize the various entries of Time-stamped Annotations
% Lists (TALs) into a structure called new_annotation

new_annotation = struct;
num_annot = 0;

for i = 1:length(annotation)
    currentstr = annotation{i};
    parselist = strsplit(currentstr, {char(43), char(20)});
    validflag = ~cellfun(@isempty, parselist);
    
    if sum(validflag) > 1 % ignore the default TALs for marking data records
        num_annot = num_annot + 1;
        
        % prep for current annotation
        firstidx = find(validflag, 1, 'first');
        onsetstr = parselist{firstidx};
        if contains(onsetstr, char(21))
            secondparse = strsplit(onsetstr, char(21));
            onsetstr = secondparse{1};
            offset = str2num(onsetstr);
            new_annotation(num_annot).offset = offset;
            new_annotation(num_annot).duration = str2num(secondparse{2});
        else
            offset = str2num(onsetstr);
            new_annotation(num_annot).offset = offset;
            new_annotation(num_annot).duration = 0;
        end
        
        % create annotation string
        contentidx = find(validflag==1);
        contentidx(contentidx == firstidx) = [];
        s = parselist(contentidx);
        annot_content = [sprintf('%s; ',s{1:end-1}),s{end}];
        
        new_annotation(num_annot).annotation = annot_content;
    end
end

end
