function [ subData ] = anabehav_P300(subID, filepath)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - P300**
%
% - Auditory Oddball Task (P300) 
%
% - loads a text file with P300 Task performance saved by a NBS
% Presentation experiment, plots basic behavioral performance, and saves a 
% single .mat file with organized behavioral performance data
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for the
%                           P300 Task txt file in the filepath directory:
%                               - subID_night2_P300.txt
%                           if the file is missing, an error will be
%                           reported and the function will terminate. 
%
%           - filepath:     a string specifying the directory containing
%                           the P300 Task txt file.
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
textfn = fullfile(filepath, [subID, '_night2_P300.txt']);
assert(isfile(textfn), [subID, '_night2_P300.txt file is missing from the directory!'])

%% Load data files
% Start loading...

% load from text file
delimiterIn = '\t';
headerlinesIn = 1;
A = importdata(textfn, delimiterIn, headerlinesIn);
data = A.data;

% prepare data into a structure
subData = struct;
sID = string;
sID(1:size(data,1),:) = subID;
subData.subID = sID;
subData.blockNum = data(:,1);
subData.trialNum = data(:,2);
subData.Condition = data(:,5); % 0 = frequent distractor, 1 = infrequent oddball target
subData.accuracy = data(:,6);
subData.RT = data(:,4);
subData.rawtxt = A;

% save the structure 
save(fullfile(filepath, [subID, '_P300_data']), 'subData')

%% Behavioral analysis

% Acc
acc_distractor = mean(subData.accuracy(subData.Condition==0));
acc_target = mean(subData.accuracy(subData.Condition==1));

% SDT measures:
FA = sum(subData.accuracy(subData.Condition==0)==0);
MISS = sum(subData.accuracy(subData.Condition==1)==0);

% Response time (RT)
rt_target = mean(subData.RT(subData.Condition==1));

% Make a plot 
figure 
set(gcf, 'Position', [100, 100, 1000, 700])
hold on
bar([rt_target])
ylabel('Response time (ms)')
yyaxis right
bar([0; 0; acc_distractor; acc_target])
ylabel('Accuracy')
xticks([1, 3, 4])
xticklabels({'Target RT', 'Distractor', 'Target'})
plot([2,2], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3)
l1=plot(NaN, NaN, '-', 'Color', 'k');
l2=plot(NaN, NaN, '-', 'Color', 'k');
legend([l1, l2], {['False Alarm: ', num2str(FA)], ['Miss: ', num2str(MISS)]}, 'Location', 'southeast')
title('P300 Task')
set(gca, 'FontSize', 20)

% save this figure
saveas(gcf, fullfile(filepath, [subID, '_P300_plot.png']))
close all

end

