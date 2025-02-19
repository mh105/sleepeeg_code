function [] = SleepEEG_vet(subID)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - VET**
%
% - used to vet all files that should be present in a subject folder and
% checks over all file naming conventions. Produce a subID_vet_info.txt
% file in the log folder that summarizes missing files
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for
%                           relevant files in the subject folder
%                           Filenames to vet are hardcoded.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - no output for this function.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Addpath to appropriate folders
% matlabroot will help detect the current environment
dataDir = SleepEEG_addpath(matlabroot);
subDir = fullfile(dataDir, subID);

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG_vet: ';

%% Open a new subID_vet_info.txt
infofn = fullfile(subDir, 'log', [subID, '_vet_info.txt']);
if isfile(infofn)
    delete(infofn)
end

diary(infofn)
disp([mHead, 'Checking subject: ', subID])
startdatetime = string(datetime);
startdatetime = strrep(startdatetime,':','_');
disp(startdatetime)
disp(' ')

%% Check all files are present for the subject
disp('[Checking all subfolders have the correct files...]')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([mHead, 'folder: ', subID, '/channel'])
p = fullfile(subDir, 'channel');
report_exist(p, [subID, '_night1_Sleep_Channel_25_all_segment.mat'])
report_exist(p, [subID, '_night1_Sleep_Channel_84_all_segment.mat'])
report_exist(p, [subID, '_night2_Sleep_Channel_25_all_segment.mat'])
report_exist(p, [subID, '_night2_Sleep_Channel_84_all_segment.mat'])
disp('completed.')

disp('--------------------------------------------------')
disp([mHead, 'folder: ', subID, '/clinical'])
p = fullfile(subDir, 'clinical');
report_exist(p, [subID, '_night1_aligned.edf'])
report_exist(p, [subID, '_night1_rvsalign.mat'])
report_exist(p, [subID, '_night1_Sleep_clinical_deidentified.edf'])
report_exist(p, [subID, '_night2_aligned.edf'])
report_exist(p, [subID, '_night2_rvsalign.mat'])
report_exist(p, [subID, '_night2_Sleep_clinical_deidentified.edf'])
disp('completed.')

disp('--------------------------------------------------')
disp([mHead, 'folder: ', subID, '/fastscan'])
p = fullfile(subDir, 'fastscan');
report_exist(p, [subID, '_night1_fastscan.mat'])
report_exist(p, [subID, '_night1_fastscan.txt'])
report_exist(p, [subID, '_night1_labelled.fsn'])
report_exist(p, [subID, '_night1.fsn'])
report_exist(p, [subID, '_night2_fastscan.mat'])
report_exist(p, [subID, '_night2_fastscan.txt'])
report_exist(p, [subID, '_night2_labelled.fsn'])
report_exist(p, [subID, '_night2.fsn'])
disp('completed.')

disp('--------------------------------------------------')
disp([mHead, 'folder: ', subID, '/raw'])
p = fullfile(subDir, 'raw');
report_exist(p, [subID, '_night1_PVT.cnt'])
report_exist(p, [subID, '_night1_PVT.evt'])
report_exist(p, [subID, '_night1_Resting.cnt'])
report_exist(p, [subID, '_night1_Resting.evt'])
report_exist(p, [subID, '_night1_ObjLoc.cnt'])
report_exist(p, [subID, '_night1_ObjLoc.evt'])
report_exist(p, [subID, '_night1_PPT_mock.cnt'])
report_exist(p, [subID, '_night1_PPT_mock.evt'])
report_exist(p, [subID, '_night1_Eyeclose.cnt'])
report_exist(p, [subID, '_night1_Eyeclose.evt'])
report_exist(p, [subID, '_night1_Sleep.cnt'])
report_exist(p, [subID, '_night1_Sleep.evt'])
report_exist(p, [subID, '_morning1_PVT.cnt'])
report_exist(p, [subID, '_morning1_PVT.evt'])
report_exist(p, [subID, '_morning1_PPT_mock.cnt'])
report_exist(p, [subID, '_morning1_PPT_mock.evt'])
report_exist(p, [subID, '_night2_PVT.cnt'])
report_exist(p, [subID, '_night2_PVT.evt'])
report_exist(p, [subID, '_night2_Resting.cnt'])
report_exist(p, [subID, '_night2_Resting.evt'])
report_exist(p, [subID, '_night2_P300.cnt'])
report_exist(p, [subID, '_night2_P300.evt'])
report_exist(p, [subID, '_night2_PPT.cnt'])
report_exist(p, [subID, '_night2_PPT.evt'])
report_exist(p, [subID, '_night2_Eyeclose.cnt'])
report_exist(p, [subID, '_night2_Eyeclose.evt'])
report_exist(p, [subID, '_night2_Sleep.cnt'])
report_exist(p, [subID, '_night2_Sleep.evt'])
report_exist(p, [subID, '_morning2_PVT.cnt'])
report_exist(p, [subID, '_morning2_PVT.evt'])
report_exist(p, [subID, '_morning2_PPT.cnt'])
report_exist(p, [subID, '_morning2_PPT.evt'])
if isfile(fullfile(p, [subID, '_night1_MST_V1.cnt']))
    report_exist(p, [subID, '_night1_MST_V1.cnt'])
    report_exist(p, [subID, '_night1_MST_V1.evt'])
    report_exist(p, [subID, '_morning1_MST_V1.cnt'])
    report_exist(p, [subID, '_morning1_MST_V1.evt'])
    report_exist(p, [subID, '_night2_MST_V2.cnt'])
    report_exist(p, [subID, '_night2_MST_V2.evt'])
    report_exist(p, [subID, '_morning2_MST_V2.cnt'])
    report_exist(p, [subID, '_morning2_MST_V2.evt'])
elseif isfile(fullfile(p, [subID, '_night1_MST_V2.cnt']))
    report_exist(p, [subID, '_night1_MST_V2.cnt'])
    report_exist(p, [subID, '_night1_MST_V2.evt'])
    report_exist(p, [subID, '_morning1_MST_V2.cnt'])
    report_exist(p, [subID, '_morning1_MST_V2.evt'])
    report_exist(p, [subID, '_night2_MST_V1.cnt'])
    report_exist(p, [subID, '_night2_MST_V1.evt'])
    report_exist(p, [subID, '_morning2_MST_V1.cnt'])
    report_exist(p, [subID, '_morning2_MST_V1.evt'])
else
    disp('ERROR finding MST files! Please manually check!')
end
disp('completed.')

disp('--------------------------------------------------')
disp([mHead, 'folder: ', subID, '/set'])
p = fullfile(subDir, 'set');
report_exist(p, [subID, '_night1_Resting_ds500_Z3.set'])
report_exist(p, [subID, '_night1_Eyeclose_ds500_Z3.set'])
report_exist(p, [subID, '_night1_Sleep_ds500_Z3.set'])
report_exist(p, [subID, '_night2_Resting_ds500_Z3.set'])
report_exist(p, [subID, '_night2_Eyeclose_ds500_Z3.set'])
report_exist(p, [subID, '_night2_Sleep_ds500_Z3.set'])
disp('completed.')

disp('--------------------------------------------------')
disp([mHead, 'folder: ', subID, '/task'])
p = fullfile(subDir, 'task');
report_exist(p, [subID, '_visit1_PVT.txt'])
report_exist(p, [subID, '_night1_PVT.txt'])
report_exist(p, [subID, '_morning1_PVT.txt'])
report_exist(p, [subID, '_night2_PVT.txt'])
report_exist(p, [subID, '_morning2_PVT.txt'])

report_exist(p, [subID, '_night2_P300.txt'])

report_exist(p, [subID, '_night1_ObjLoc.txt'])

report_exist(p, [subID, '_night1_2pairs_Picture_Pair_mock_Results_Test1.txt'])
report_exist(p, [subID, '_night1_2pairs_Picture_Pair_mock_StudyPhase.txt'])
report_exist(p, [subID, '_night1_2pairs_Picture_Pair_mock_Test2_List.txt'])
report_exist(p, [subID, '_night1_4pairs_Picture_Pair_mock_Results_Test1.txt'])
report_exist(p, [subID, '_night1_4pairs_Picture_Pair_mock_Results_Test2.txt'])
report_exist(p, [subID, '_night1_4pairs_Picture_Pair_mock_StudyPhase.txt'])
report_exist(p, [subID, '_night1_4pairs_Picture_Pair_mock_Test2_List.txt'])
report_exist(p, [subID, '_night2_Picture_Pair_Results_Test1.txt'])
report_exist(p, [subID, '_night2_Picture_Pair_Results_Test2.txt'])
report_exist(p, [subID, '_night2_Picture_Pair_StudyPhase.txt'])
report_exist(p, [subID, '_night2_Picture_Pair_Test2_List.txt'])

if isfile(fullfile(p, [subID, '_night1_MST_V1_Test1.txt']))
    report_exist(p, [subID, '_night1_MST_V1_Test1.txt'])
    report_exist(p, [subID, '_morning1_MST_V1_Test2.txt'])
    report_exist(p, [subID, '_night2_MST_V2_Test1.txt'])
    report_exist(p, [subID, '_morning2_MST_V2_Test2.txt'])
elseif isfile(fullfile(p, [subID, '_night1_MST_V2_Test1.cnt']))
    report_exist(p, [subID, '_night1_MST_V2_Test1.txt'])
    report_exist(p, [subID, '_morning1_MST_V2_Test2.txt'])
    report_exist(p, [subID, '_night2_MST_V1_Test1.txt'])
    report_exist(p, [subID, '_morning2_MST_V1_Test2.txt'])
else
    disp('ERROR finding MST files! Please manually check!')
end
disp('completed.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Check for signal quality of EEG .cnt files
disp(' ')
disp(' ')
disp('[Checking the number of triggers and multiple segments in raw .cnt files...]')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp([mHead, 'folder: ', subID, '/raw'])
p = fullfile(subDir, 'raw');
disp(' ')
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_PVT.cnt']);
% this is administered on iPad, no trigger
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,PVT'), {EEG1.event.type}))
        disp('Cannot find "1000,PVT" task annotation.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_Resting.cnt']);
% there should be triggers marking start/end of eyes open and eyes closed
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,Resting'), {EEG1.event.type}))
        disp('Cannot find "1000,Resting" task annotation.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1001,EyesOpen_Start'), {EEG1.event.type}))
        disp('Cannot find "EyesOpen_Start" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1002,EyesOpen_End'), {EEG1.event.type}))
        disp('Cannot find "EyesOpen_End" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1003,EyesClose_Start'), {EEG1.event.type}))
        disp('Cannot find "EyesClose_Start" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1004,EyesClose_End'), {EEG1.event.type}))
        disp('Cannot find "EyesClose_End" trigger.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_ObjLoc.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,ObjLoc'), {EEG1.event.type}))
        disp('Cannot find "1000,ObjLoc" task annotation.')
    end
    obj_num = 0;
    loc_num = 0;
    both_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '3')
            obj_num = obj_num + 1;
        elseif strcmp(EEG1.event(ii).type, '4')
            loc_num = loc_num + 1;
        elseif strcmp(EEG1.event(ii).type, '5')
            both_num = both_num + 1;
        end
    end
    if obj_num ~= 24 || loc_num ~= 24 || both_num ~= 24
        disp(['There are triggers for ', num2str(obj_num), ' object trials, ',...
            num2str(loc_num), ' location trials, and ',...
            num2str(both_num), ' both trials in the Object Location Task.'])
    else
        disp('Successfully found 24 trials for each condition in the Object Location Task.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_PPT_mock.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'mock', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "PPT mock" task annotation.')
    end
    study_num = 0;
    test_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '2')
            study_num = study_num + 1;
        elseif strcmp(EEG1.event(ii).type, '3')
            test_num = test_num + 1;
        end
    end
    disp(['There are triggers for ', num2str(study_num), ' study pairs, and ',...
        num2str(test_num), ' immediate test pairs in the Picture Pairing mock Task.'])
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_Eyeclose.cnt']);
% no trigger for this task
if ~isempty(EEG1)
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night1_Sleep.cnt']);
% variable number of triggers, just check for segment.
if ~isempty(EEG1)
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_morning1_PVT.cnt']);
% this is administered on iPad, no trigger
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,PVT'), {EEG1.event.type}))
        disp('Cannot find "1000,PVT" task annotation.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_morning1_PPT_mock.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'mock', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "PPT mock" task annotation.')
    end
    test_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '3')
            test_num = test_num + 1;
        end
    end
    if test_num ~= 12
        disp(['There are triggers for ', num2str(test_num), ' delayed test pairs in the Picture Pairing mock Task.'])
    else
        disp('Successfully found 12 delayed test pairs in the Picture Pairing mock Task.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_PVT.cnt']);
% this is administered on iPad, no trigger
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,PVT'), {EEG1.event.type}))
        disp('Cannot find "1000,PVT" task annotation.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_Resting.cnt']);
% there should be triggers marking start/end of eyes open and eyes closed
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,Resting'), {EEG1.event.type}))
        disp('Cannot find "1000,Resting" task annotation.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1001,EyesOpen_Start'), {EEG1.event.type}))
        disp('Cannot find "EyesOpen_Start" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1002,EyesOpen_End'), {EEG1.event.type}))
        disp('Cannot find "EyesOpen_End" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1003,EyesClose_Start'), {EEG1.event.type}))
        disp('Cannot find "EyesClose_Start" trigger.')
    end
    if ~any(cellfun(@(x) strcmp(strrep(x,' ',''), '1004,EyesClose_End'), {EEG1.event.type}))
        disp('Cannot find "EyesClose_End" trigger.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_P300.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,P300'), {EEG1.event.type}))
        disp('Cannot find "1000,P300" task annotation.')
    end
    odd_num = 0;
    reg_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '2')
            reg_num = reg_num + 1;
        elseif strcmp(EEG1.event(ii).type, '3')
            odd_num = odd_num + 1;
        end
    end
    if reg_num ~= 160 || odd_num ~= 40
        disp(['There are triggers for ', num2str(reg_num), ' distractor tones, and ',...
            num2str(odd_num), ' target tones in the P300 Task.'])
    else
        disp('Successfully found 160 distractor tones and 40 target tones in the P300 Task.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_PPT.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,PPT'), {EEG1.event.type}))
        disp('Cannot find "1000,PPT" task annotation.')
    end
    study_num = 0;
    test_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '2')
            study_num = study_num + 1;
        elseif strcmp(EEG1.event(ii).type, '3')
            test_num = test_num + 1;
        end
    end
    if study_num ~= 80 || test_num ~= 60
        disp(['There are triggers for ', num2str(study_num), ' study pairs, and ',...
            num2str(test_num), ' immediate test pairs in the Picture Pairing Task.'])
    else
        disp('Successfully found 80 study pairs and 60 immediate test pairs in the Picture Pairing Task.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_Eyeclose.cnt']);
% no trigger for this task
if ~isempty(EEG1)
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_night2_Sleep.cnt']);
% variable number of triggers, just check for segment.
if ~isempty(EEG1)
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_morning2_PVT.cnt']);
% this is administered on iPad, no trigger
if ~isempty(EEG1)
    if ~any(cellfun(@(x) strcmpi(strrep(x,' ',''), '1000,PVT'), {EEG1.event.type}))
        disp('Cannot find "1000,PVT" task annotation.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, [subID, '_morning2_PPT.cnt']);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'PPT', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "1000,PPT" task annotation.')
    end
    test_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '3')
            test_num = test_num + 1;
        end
    end
    if test_num ~= 60
        disp(['There are triggers for ', num2str(test_num), ' delayed test pairs in the Picture Pairing Task.'])
    else
        disp('Successfully found 60 delayed test pairs in the Picture Pairing Task.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
if isfile(fullfile(p, [subID, '_night1_MST_V1.cnt']))
    v1t1 = [subID, '_night1_MST_V1.cnt'];
    v1t2 = [subID, '_morning1_MST_V1.cnt'];
    v2t1 = [subID, '_night2_MST_V2.cnt'];
    v2t2 = [subID, '_morning2_MST_V2.cnt'];
elseif isfile(fullfile(p, [subID, '_night1_MST_V2.cnt']))
    v1t1 = [subID, '_night2_MST_V1.cnt'];
    v1t2 = [subID, '_morning2_MST_V1.cnt'];
    v2t1 = [subID, '_night1_MST_V2.cnt'];
    v2t2 = [subID, '_morning1_MST_V2.cnt'];
else
    disp('ERROR finding MST files! Please manually check!')
    disp(' ')
    v1t1 = [subID, '_night1_MST_V1.cnt'];
    v1t2 = [subID, '_morning1_MST_V1.cnt'];
    v2t1 = [subID, '_night2_MST_V2.cnt'];
    v2t2 = [subID, '_morning2_MST_V2.cnt'];
end
%%
EEG1 = check_segment(p, v1t1);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'MST', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "1000,MST" task annotation.')
    end
    start_num = 0;
    end_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '22')
            start_num = start_num + 1;
        elseif strcmp(EEG1.event(ii).type, '23')
            end_num = end_num + 1;
        end
    end
    if start_num ~= 12 || end_num~= 12
        disp(['There are triggers for ', num2str(start_num), ' start and ',...
            num2str(end_num), ' end of blocks in the MST V1_Test1.'])
    else
        disp('Successfully found 12 blocks in the MST V1_Test1.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, v1t2);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'MST', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "1000,MST" task annotation.')
    end
    start_num = 0;
    end_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '22')
            start_num = start_num + 1;
        elseif strcmp(EEG1.event(ii).type, '23')
            end_num = end_num + 1;
        end
    end
    if start_num ~= 12 || end_num~= 12
        disp(['There are triggers for ', num2str(start_num), ' start and ',...
            num2str(end_num), ' end of blocks in the MST V1_Test2.'])
    else
        disp('Successfully found 12 blocks in the MST V1_Test2.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, v2t1);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'MST', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "1000,MST" task annotation.')
    end
    start_num = 0;
    end_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '22')
            start_num = start_num + 1;
        elseif strcmp(EEG1.event(ii).type, '23')
            end_num = end_num + 1;
        end
    end
    if start_num ~= 12 || end_num~= 12
        disp(['There are triggers for ', num2str(start_num), ' start and ',...
            num2str(end_num), ' end of blocks in the MST V2_Test1.'])
    else
        disp('Successfully found 12 blocks in the MST V2_Test1.')
    end
    disp(' ')
    disp('completed.')
end
disp('--------------------------------------------------')
%%
EEG1 = check_segment(p, v2t2);
if ~isempty(EEG1)
    if ~any(cellfun(@(x) contains(strrep(x,' ',''), 'MST', 'IgnoreCase',true), {EEG1.event.type}))
        disp('Cannot find "1000,MST" task annotation.')
    end
    start_num = 0;
    end_num = 0;
    for ii = 1:length(EEG1.event)
        if strcmp(EEG1.event(ii).type, '22')
            start_num = start_num + 1;
        elseif strcmp(EEG1.event(ii).type, '23')
            end_num = end_num + 1;
        end
    end
    if start_num ~= 12 || end_num~= 12
        disp(['There are triggers for ', num2str(start_num), ' start and ',...
            num2str(end_num), ' end of blocks in the MST V2_Test2.'])
    else
        disp('Successfully found 12 blocks in the MST V2_Test2.')
    end
    disp(' ')
    disp('completed.')
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%
diary off
close all

end

%% HELPER FUNCTIONS
function [] = report_exist(fpath, fname)
if ~isfile(fullfile(fpath, fname))
    disp(['Missing file: [ ', fname, ' ]'])
end
end

function [ EEG1 ] = check_segment(fpath, fname)
disp(' ')
disp(['Checking the .cnt file ', fname, '...'])
disp(' ')

if isfile(fullfile(fpath, fname))
    diary off
    EEG1 = ANT_interface_readcnt(fname, fpath, [true, 500], true);
    diary on
    % Check impedance
    if ~strcmp(EEG1.event(1).type, '0, Impedance')
        disp('Start impedance measures are missing!')
    end
    if ~strcmp(EEG1.event(end).type, '0, Impedance')
        disp('End impedance measures are missing!')
    end
    % Check number of segments and interruptions
    disconnect_num = 0;
    add_imp_num = 0;
    for ii = 2:length(EEG1.event)-1
        if strcmp(EEG1.event(ii).type, '0, Impedance')
            add_imp_num = add_imp_num + 1;
        elseif strcmp(EEG1.event(ii).type, '9001, Amplifier disconnected')
            disconnect_num = disconnect_num + 1;
        end
    end
    if disconnect_num > 0 || add_imp_num > 0
        disp([num2str(add_imp_num), ' additional impedance checks found during recording.'])
        disp(['Number of Amplifier disconnection message found: ', num2str(disconnect_num)])
        disp(['Recording is broken into ', num2str(disconnect_num+1), ' segments.'])
    end
else
    disp(['Missing file: [ ', fname, ' ]'])
    EEG1 = [];
end
end
