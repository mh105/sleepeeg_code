function [ subData ] = anabehav_PDP(subID, filepath)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - PDP**
%
% - Process Dissociation Paradigm (PDP)
%
% - Also called Picture Pairing Task or Episodic Memory Task
%
% - loads two text files with Process Dissociation Paradigm performance
% saved by NBS Presentation experiments, plots basic behavioral
% performance, and saves a single .mat file with organized behavioral
% performance data
%
% - Edit note: need to add a section to do sanity check on the mock version
% performance. Right now we don't have the data so we will skip it.
%
% **Need to add codes related to the mock test! 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for the
%                           PDP txt files in the filepath directory:
%                               - subID_night2_Picture_Pair_Results_V1_Test1.txt
%                               - subID_night2_Picture_Pair_Results_V1_Test2.txt
%                           if any of the files is missing, an error will be
%                           reported and the function will terminate.
%                               - mock version filenames are variable, and
%                               will be searched within the same directory
%
%           - filepath:     a string specifying the directory containing
%                           the Process Dissociation Paradigm txt files.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - subData:      a structure containing task variables and
%                           performances for the analyzed subject. This
%                           variable is also saved to the filepath as a
%                           .mat file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% file validity check
try
    test1fn = fullfile(filepath, [subID, '_night2_Picture_Pair_Results_V1_Test1.txt']);
    test2fn = fullfile(filepath, [subID, '_night2_Picture_Pair_Results_V1_Test2.txt']);
    assert(isfile(test1fn), [subID, '_night2_Picture_Pair_Results_V1_Test1.txt file is missing from the directory!'])
    assert(isfile(test2fn), [subID, '_night2_Picture_Pair_Results_V1_Test2.txt file is missing from the directory!'])
catch
    test1fn = fullfile(filepath, [subID, '_night2_Picture_Pair_Results_Test1.txt']);
    test2fn = fullfile(filepath, [subID, '_night2_Picture_Pair_Results_Test2.txt']);
    assert(isfile(test1fn), [subID, '_night2_Picture_Pair_Results_Test1.txt file is missing from the directory!'])
    assert(isfile(test2fn), [subID, '_night2_Picture_Pair_Results_Test2.txt file is missing from the directory!'])
end

% sanity check
testcheck = importdata(test1fn);
assert(strcmp(testcheck.textdata{1,1}, 'Test1'), 'Something is wrong with the exported _Test1.txt file.')
testcheck = importdata(test2fn);
assert(strcmp(testcheck.textdata{1,1}, 'Test2'), 'Something is wrong with the exported _Test2.txt file.')

%% Load data files
% Start loading...

opts = delimitedTextImportOptions("NumVariables", 10);
% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Trial", "Pic", "Resp", "RT", "Type", "corr", "LeftResp", "LeftRT", "RightResp", "RightRT"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
PicturePairResultsV1Test1 = readtable(test1fn, opts);
PicturePairResultsV1Test2 = readtable(test2fn, opts);
% Clear temporary variables
clear opts

% sanity check
assert(size(PicturePairResultsV1Test1,1) == 60 && size(PicturePairResultsV1Test1,2) == 10, 'Something is wrong with the dimension of imported _Test1.txt file.')
assert(size(PicturePairResultsV1Test2,1) == 60 && size(PicturePairResultsV1Test2,2) == 10, 'Something is wrong with the dimension of imported _Test2.txt file.')

% prepare data into a structure
subData = struct;
sID = string;
sID(1:size(PicturePairResultsV1Test1,1)*2,:) = subID;
subData.subID = sID;
subData.delayed = [zeros(size(PicturePairResultsV1Test1,1), 1); ones(size(PicturePairResultsV1Test2,1), 1)];
subData.trialNum = [PicturePairResultsV1Test1.Trial; PicturePairResultsV1Test2.Trial];
subData.picID = [PicturePairResultsV1Test1.Pic; PicturePairResultsV1Test2.Pic];
subData.Condition = [PicturePairResultsV1Test1.Type; PicturePairResultsV1Test2.Type]; % 0 = intact, 1 = rearranged, 2 = novel
subData.Response = [PicturePairResultsV1Test1.Resp; PicturePairResultsV1Test2.Resp];
subData.accuracy = [PicturePairResultsV1Test1.corr; PicturePairResultsV1Test2.corr];
subData.RT = [PicturePairResultsV1Test1.RT; PicturePairResultsV1Test2.RT];
rawtxt = struct;
rawtxt.immediate_test = PicturePairResultsV1Test1;
rawtxt.delayed_test = PicturePairResultsV1Test2;
subData.rawtxt = rawtxt;

% save the structure
save(fullfile(filepath, [subID, '_PDP_data']), 'subData')

%% Behavioral analysis

% Plot the performance
t1intact = sum(PicturePairResultsV1Test1.corr(PicturePairResultsV1Test1.Type==0));
t1rearrg = sum(PicturePairResultsV1Test1.corr(PicturePairResultsV1Test1.Type==1));
t1novel  = sum(PicturePairResultsV1Test1.corr(PicturePairResultsV1Test1.Type==2));

t2intact = sum(PicturePairResultsV1Test2.corr(PicturePairResultsV1Test2.Type==0));
t2rearrg = sum(PicturePairResultsV1Test2.corr(PicturePairResultsV1Test2.Type==1));
t2novel  = sum(PicturePairResultsV1Test2.corr(PicturePairResultsV1Test2.Type==2));

% some sanity checks 
assert(t1intact == sum(subData.accuracy == 1 & subData.Condition == 0 & subData.delayed == 0), 'Sanity check failed!')
assert(t1rearrg == sum(subData.accuracy == 1 & subData.Condition == 1 & subData.delayed == 0), 'Sanity check failed!')
assert(t1novel == sum(subData.accuracy == 1 & subData.Condition == 2 & subData.delayed == 0), 'Sanity check failed!')
assert(t2intact == sum(subData.accuracy == 1 & subData.Condition == 0 & subData.delayed == 1), 'Sanity check failed!')
assert(t2rearrg == sum(subData.accuracy == 1 & subData.Condition == 1 & subData.delayed == 1), 'Sanity check failed!')
assert(t2novel == sum(subData.accuracy == 1 & subData.Condition == 2 & subData.delayed == 1), 'Sanity check failed!')

% behavioral performance plot 
figure
ax = figdesign(2,2,'merge',{3:4});
for ii = 1:length(ax); title(ax(ii), ii); end
set(gcf,'units', 'pixels', 'Position', [0, 0, 1100, 1000])
axes(ax(1));
hold on
bar([t1intact, t1rearrg, t1novel])
plot([0,4], [10,10], 'r--', 'LineWidth', 3)
yticks(1:20)
xticks([1, 2, 3])
xticklabels({'Intact', 'Rearrange', 'Novel'})
title('Immediate Test')
xlabel('Trial Conditions')
ylabel('# Correct Responses')
ylim([0, 20])
set(gca, 'FontSize', 20)

axes(ax(2));
hold on
bar([t2intact, t2rearrg, t2novel])
plot([0,4], [10,10], 'r--', 'LineWidth', 3)
yticks(1:20)
xticks([1, 2, 3])
xticklabels({'Intact', 'Rearrange', 'Novel'})
title('Delayed Test')
xlabel('Trial Conditions')
ylabel('# Correct Responses')
ylim([0, 20])
set(gca, 'FontSize', 20)

subjectID = strrep(subID,'_','-');
annotstring = ['SubID: ', subjectID];
annotation('textbox',[0.475, 0.79, 0.1, 0.2],'String',annotstring,'FitBoxToText','on')

% Behavioral in-depth visualization 
% visualize accuracy and RT on the same plot 

trialN = [1:60, 66:125];
axes(ax(3));
hold on
c = plot(trialN(subData.accuracy == 1), subData.RT(subData.accuracy == 1)/1000, '.', 'Color', [0, 0.75, 0], 'MarkerSize', 20);
ic_intact = plot(trialN(subData.accuracy == 0 & subData.Condition == 0), subData.RT(subData.accuracy == 0 & subData.Condition == 0)/1000, '*', 'Color', [0, 0, 0.75], 'MarkerSize', 20);
ic_rearrange = plot(trialN(subData.accuracy == 0 & subData.Condition == 1), subData.RT(subData.accuracy == 0 & subData.Condition == 1)/1000, 'h', 'Color', [0.75, 0, 0], 'MarkerSize', 15, 'MarkerFaceColor', [0.75, 0, 0]);
ic_novel = plot(trialN(subData.accuracy == 0 & subData.Condition == 2), subData.RT(subData.accuracy == 0 & subData.Condition == 2)/1000, 'x', 'Color', [0.75, 0, 0], 'MarkerSize', 20);
plot([63,63], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)
xticks([30, 96])
xticklabels({'Immediate Test 60 trials', 'Delayed Test 60 trials'})
ylabel('Response Time (s)')
title('PDP RT and Accuracy')
legend([c, ic_intact, ic_rearrange, ic_novel], {'Correct', 'Inc-Intact', 'Inc-Rearrg', 'Inc-Novel'}, 'Location', 'best')
set(gca, 'FontSize', 20)
ylimit = ylim;
axis tight 
ylim(ylimit)

% save this figure
saveas(gcf, fullfile(filepath, [subID, '_PDP_plot.png']))
close all

%% %% Compute recollection/familiarity estimate under Process Dissociation Process

% This section is commented out for now, until I discuss with Jess and Brad
% about the correct method to compute the recollection/familiarity
% estimates. 

% % store values in table
% PP_task.subID{ii} = subjectID;
% 
% PP_task.Imm_intact(ii) = t1intact;
% PP_task.Imm_rearrg(ii) = t1rearrg;
% PP_task.Imm_novel(ii) = t1novel;
% PP_task.Del_intact(ii) = t2intact;
% PP_task.Del_rearrg(ii) = t2rearrg;
% PP_task.Del_novel(ii) = t2novel;
% 
% PP_task.Imm_PI(ii) = t1intact/20;
% PP_task.Imm_PR(ii) = (20-t1rearrg)/20;
% PP_task.Imm_PN(ii) = (20-t1novel)/20;
% PP_task.Del_PI(ii) = t2intact/20;
% PP_task.Del_PR(ii) = (20-t2rearrg)/20;
% PP_task.Del_PN(ii) = (20-t2novel)/20;
% 
% %%
% % 'Wolk' method
% [R1, F1] = compute_RF(t1intact/20, (20-t1rearrg)/20, (20-t1novel)/20, 'Wolk');
% [R2, F2] = compute_RF(t2intact/20, (20-t2rearrg)/20, (20-t2novel)/20, 'Wolk');
% 
% PP_task.Imm_Wolk_R(ii) = R1;
% PP_task.Imm_Wolk_F(ii) = F1;
% PP_task.Del_Wolk_R(ii) = R2;
% PP_task.Del_Wolk_F(ii) = F2;
% 
% % 'Pair' method
% [R1, F1] = compute_RF(t1intact/20, (20-t1rearrg)/20, (20-t1novel)/20, 'Pair');
% [R2, F2] = compute_RF(t2intact/20, (20-t2rearrg)/20, (20-t2novel)/20, 'Pair');
% 
% fprintf('<strong> Night Session: </strong> \n')
% disp(['Recollection: ', num2str(R1), ' prob; Familiarity: ' num2str(F1) , ' std'])
% 
% fprintf('<strong> Morning Session: </strong> \n')
% disp(['Recollection: ', num2str(R2), ' prob; Familiarity: ' num2str(F2) , ' std'])
% 
% PP_task.Imm_Pair_R(ii) = R1;
% PP_task.Imm_Pair_F(ii) = F1;
% PP_task.Del_Pair_R(ii) = R2;
% PP_task.Del_Pair_F(ii) = F2;

end

