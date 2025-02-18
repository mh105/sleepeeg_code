function [] = SleepEEG_report(subID, project)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - REPORT**
%
% - used to report the available .cnt files under the subID directory and
% the length of each EEG recording
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 2
    project = '';
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('---------------------------')
disp([mHead 'SleepEEG_report()']);
disp('---------------------------')

%% Loop through all available .cnt files
[ ~, ~, fnsuffixlist, ~ ] = SleepEEG_configDir(subID, [], false, project);

disp(' ')
disp("%%-------------------------------%%%---------------------------------%%")
disp([mHead, 'Reporting all available .cnt files for subject <<< ' subID ' >>>:'])

for i = 1:length(fnsuffixlist)
    [dataDir, datafn, fileID, ~] = SleepEEG_configDir(subID, fnsuffixlist{i}, false, project);
    cnt_info = eepv4_read_info(datafn);
    tot_sample_point = cnt_info.sample_count;
    record_time = cnt_info.sample_count/cnt_info.sample_rate/60;

    fprintf(mSpace);
    disp(['[ ', num2str(i), ' ]  ', fileID, '.cnt'])
    disp([mSpace, 'Total recording duration = ' num2str(record_time) ' minutes, ', num2str(record_time/60) ' hours.'])
    disp([mSpace, 'Total sample points = ' num2str(tot_sample_point) '.'])
    disp([mSpace, 'Original collection sampling rate = ' num2str(cnt_info.sample_rate) 'Hz.'])

    if isfile(fullfile(dataDir, subID, 'set', [subID, '_', fnsuffixlist{i}, '_ds500_Z3', '.set']))
        disp([mSpace, 'Downsampled to 500Hz = ', subID, '_', fnsuffixlist{i}, '_ds500_Z3', '.set'])
    else
        disp([mSpace, 'Downsampled to 500Hz = not found.'])
    end

    disp(' ')
end
disp("%%-------------------------------%%%---------------------------------%%")
disp(' ')

end
