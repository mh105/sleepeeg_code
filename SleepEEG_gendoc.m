%
% **ADSLEEPEEG_PREPROCESSING SCRIPT - GENDOC**
%
% - this is the function used to generate documentations of all other
% functions. Not used for processing
%
%%
diary('Functionheader.txt')

function_list = dir('SleepEEG_*');

for ii = 1:length(function_list)
    eval(['help ', char(function_list(ii).name)])
    disp(' ')
    disp(' ')
    disp('=============================================================================')
    disp(' ')
    disp(' ')
end
diary off
