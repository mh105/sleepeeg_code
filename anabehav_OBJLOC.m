function [ subData ] = anabehav_OBJLOC(subID, filepath)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - ObjLoc**
%
% - Object Location Task (ObjLoc) 
%
% - loads a text file with Object Location Task performance saved by a NBS
% Presentation experiment, plots basic behavioral performance, and saves a
% single .mat file with organized behavioral performance data
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for the
%                           ObjLoc Task txt file in the filepath directory:
%                               - subID_night1_ObjLoc.txt
%                           if the file is missing, an error will be
%                           reported and the function will terminate. 
%
%           - filepath:     a string specifying the directory containing
%                           the Object Location Task txt file.
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
textfn = fullfile(filepath, [subID, '_night1_ObjLoc.txt']);
assert(isfile(textfn), [subID, '_night1_ObjLoc.txt file is missing from the directory!'])

%% Load data files
% Start loading...

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Block", "Trial", "Cond", "Studied", "Resp", "Corr", "RT"];
opts.VariableTypes = ["double", "double", "categorical", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, "Cond", "EmptyFieldRule", "auto");
% Import the data
data = readtable(textfn, opts);
% Clear temporary variables
clear opts

% prepare data into a structure
subData = struct;
sID = string;
sID(1:size(data,1),:) = subID;
subData.subID = sID;
subData.blockNum = data.Block;
subData.trialNum = data.Trial;
subData.Condition = data.Cond;
subData.Response = data.Resp;
subData.accuracy = data.Corr;
subData.RT = data.RT;
subData.rawtxt = data;

% save the structure 
save(fullfile(filepath, [subID, '_ObjLoc_data']), 'subData')

%% Behavioral analysis

% Acc
acc_obj = mean(data.Corr(data.Cond == 'Obj'));
acc_loc = mean(data.Corr(data.Cond == 'Loc'));
acc_both = mean(data.Corr(data.Cond == 'Both'));

% Response time (RT)
rt_obj = mean(data.RT(data.Cond == 'Obj'));
rt_loc = mean(data.RT(data.Cond == 'Loc'));
rt_both = mean(data.RT(data.Cond == 'Both'));

% Make a plot 
figure 
hold on
bar([rt_obj, 0; rt_loc, 0; rt_both, 0])
ylabel('Response time (ms)')
yyaxis right
bar([0 acc_obj; 0, acc_loc; 0, acc_both])
ylabel('Accuracy')
xticklabels({'Object', ' ', 'Location', ' ', 'Both'})
title('Object Location Task')
set(gca, 'FontSize', 20)

% save this figure
saveas(gcf, fullfile(filepath, [subID, '_ObjLoc_plot.png']))
close all

end

