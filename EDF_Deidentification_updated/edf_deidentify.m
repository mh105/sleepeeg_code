function edf_deidentify(filename,savename)
% Function to deidentify any patient information in the header in the EDF file.

% Input: filename in text format

% Output: deidentified file saved as savename

% Usage: Example : edf_deidentify('XYZ.edf','XYZ_deidentified.edf');

%*********************************************************************************

% load data file
[ edf_header ]=SleepEEG_loadedf(filename);

% Extract subject ID
strlist = strsplit(filename, '_');

% Extract night session number
nightidx = strfind(filename, 'night');
assert(~isempty(nightidx), 'No "night" string contained in fnsuffix!')
sessionid = ['night', num2str(filename(nightidx+5))];

% Update patient name field with subject ID and night session
edf_header.patient_id = ['X X X ', strlist{1}, '_', sessionid];

% close file explicitly
fclose all;

% write to a new edf file with _identified suffix
status = blockEdfWrite(savename, edf_header);

end