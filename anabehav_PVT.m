function [ subData ] = anabehav_PVT(subID, filepath)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - PVT**
%
% - Psychomotor Vigilance Task (PVT)
%
% - loads a set of 5 text files with PVT performance saved by NBS
% Presentation experiments administered on an iPad, plots basic behavioral
% performance, and saves a single .mat file with organized behavioral
% performance data 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for the
%                           5 PVT task files in the filepath directory:
%                               - subID_visit1_PVT.txt
%                               - subID_morning2_PVT.txt
%                               - subID_morning1_PVT.txt
%                               - subID_night1_PVT.txt
%                               - subID_night2_PVT.txt
%                           if any of these files is missing, an error will
%                           be reported and the function will terminate.
%
%           - filepath:     a string specifying the directory containing
%                           the MST task files.
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
visit1_fn = [subID, '_visit1_PVT.txt'];
night1_fn = [subID, '_night1_PVT.txt'];
morning1_fn = [subID, '_morning1_PVT.txt'];
night2_fn = [subID, '_night2_PVT.txt'];
morning2_fn = [subID, '_morning2_PVT.txt'];

assert(isfile(fullfile(filepath, visit1_fn)), [visit1_fn, ' file is missing from the directory!'])
assert(isfile(fullfile(filepath, night1_fn)), [night1_fn, ' file is missing from the directory!'])
assert(isfile(fullfile(filepath, morning1_fn)), [morning1_fn, ' file is missing from the directory!'])
assert(isfile(fullfile(filepath, night2_fn)), [night2_fn, ' file is missing from the directory!'])
assert(isfile(fullfile(filepath, morning2_fn)), [morning2_fn, ' file is missing from the directory!'])

visit1_fn = fullfile(filepath, visit1_fn);
night1_fn = fullfile(filepath, night1_fn);
morning1_fn = fullfile(filepath, morning1_fn);
night2_fn = fullfile(filepath, night2_fn);
morning2_fn = fullfile(filepath, morning2_fn);

%% Load data files
% Start loading...

% due to the report of false alarms and misses at the end of the text file,
% we can't use importdata to read the .txt file directly. We need to read
% it in as raw strings and parse/convert data types manually. 

[ v1_data, v1_fa, v1_miss, v1_raw ] = readPVTtxt(visit1_fn);
[ n1_data, n1_fa, n1_miss, n1_raw ] = readPVTtxt(night1_fn);
[ m1_data, m1_fa, m1_miss, m1_raw ] = readPVTtxt(morning1_fn);
[ n2_data, n2_fa, n2_miss, n2_raw ] = readPVTtxt(night2_fn);
[ m2_data, m2_fa, m2_miss, m2_raw ] = readPVTtxt(morning2_fn);

assert(size(v1_data,1) > 50, 'Not enough trials in _visit1_PVT.txt. Check the file and re-decrypt!')
assert(size(n1_data,1) > 50, 'Not enough trials in _night1_PVT.txt. Check the file and re-decrypt!')
assert(size(m1_data,1) > 50, 'Not enough trials in _morning1_PVT.txt. Check the file and re-decrypt!')
assert(size(n2_data,1) > 50, 'Not enough trials in _night2_PVT.txt. Check the file and re-decrypt!')
assert(size(m2_data,1) > 50, 'Not enough trials in _morning2_PVT.txt. Check the file and re-decrypt!')

totallength = size(v1_data,1) + size(n1_data,1) + size(m1_data,1) + ...
    size(n2_data,1) + size(m2_data, 1);

% Organizing into a structure variable
subData = struct;
sID = string;
sID(1:totallength,:) = subID;
subData.subID = sID;
subData.visit1 = [ones(size(v1_data,1),1); zeros(size(n1_data,1) + size(m1_data,1) + ...
    size(n2_data,1) + size(m2_data, 1), 1)];
subData.night = [zeros(size(v1_data,1),1); ones(size(n1_data,1),1)*1; zeros(size(m1_data,1),1);...
    ones(size(n1_data,1),1)*2; zeros(size(m1_data,1),1)];
subData.morning = [zeros(size(v1_data,1),1); zeros(size(n1_data,1),1); ones(size(m1_data,1),1)*1;...
    zeros(size(n1_data,1),1); ones(size(m1_data,1),1)*2];
subData.trialNum = [v1_data(:,1); n1_data(:,1); m1_data(:,1); n2_data(:,1); m2_data(:,1)];
subData.RT = [v1_data(:,2); n1_data(:,2); m1_data(:,2); n2_data(:,2); m2_data(:,2)];
subData.ISI = [v1_data(:,3); n1_data(:,3); m1_data(:,3); n2_data(:,3); m2_data(:,3)];

% store the false alarm and miss values 
SDT_count.session = {'visit1', 'night1', 'morning1', 'night2', 'morning2'};
SDT_count.false_alarms = [v1_fa, n1_fa, m1_fa, n2_fa, m2_fa];
SDT_count. misses = [v1_miss, n1_miss, m1_miss, n2_miss, m2_miss];
subData.SDT_count = SDT_count;

% also store the raw text file data 
rawtxt = struct;
rawtxt.visit1 = v1_raw;
rawtxt.night1 = n1_raw;
rawtxt.morning1 = m1_raw;
rawtxt.night2 = n2_raw;
rawtxt.morning2 = m2_raw;
subData.rawtxt = rawtxt;

% save the structure 
save(fullfile(filepath, [subID, '_PVT_data']), 'subData')

%% Behavioral analysis

% RT measures

%Plot performance over 3min across all sessions
figure
set(gcf,'Position',[500 500 1000 800])
% averaging across the two versions
subplot(2,1,1)
hold on
plot(v1_data(:,2), '.-', 'MarkerSize',15)
plot(n1_data(:,2), '.-', 'MarkerSize',15)
plot(m1_data(:,2), '.-', 'MarkerSize',15)
plot(n2_data(:,2), '.-', 'MarkerSize',15)
plot(m2_data(:,2), '.-', 'MarkerSize',15)
title('PVT Response Time across Sessions')
ylabel('Response time (ms)')
xlabel('# stimulus in 3-min blocks')
legend('Visit1', 'Night1', 'Morning1', 'Night2', 'Morning2', 'Location', 'best')
set(gca, 'FontSize', 20)

ax = subplot(2,1,2);
SDT_count.false_alarms(SDT_count.false_alarms==0) = 0.1;
SDT_count.misses(SDT_count.misses==0) = 0.1;
hold on
ymax = max([SDT_count.false_alarms, SDT_count.misses])+1;
bar([SDT_count.false_alarms(1), 0, SDT_count.false_alarms(2), 0, SDT_count.false_alarms(3), 0, ...
    SDT_count.false_alarms(4), 0, SDT_count.false_alarms(5), 0])
ylim([0, ymax])
ylabel('# False Alarms')
yyaxis right
bar([0, SDT_count.misses(1), 0, SDT_count.misses(2), 0, SDT_count.misses(3),...
    0, SDT_count.misses(4), 0, SDT_count.misses(5)])
ylim([0, ymax])
ylabel('# Misses')
xticks([1.5, 3.5, 5.5, 7.5, 9.5])
xticklabels({'Visit1', 'Night1', 'Morning1', 'Night2', 'Morning2'})
title('PVT False Alarms and Misses')
set(gca, 'FontSize', 20)
plot([2.5,2.5], ylim, '--', 'Color', [0.5667,0.5667,0.5667], 'LineWidth', 2)
plot([4.5,4.5], ylim, '--', 'Color', [0.5667,0.5667,0.5667], 'LineWidth', 2)
plot([6.5,6.5], ylim, '--', 'Color', [0.5667,0.5667,0.5667], 'LineWidth', 2)
plot([8.5,8.5], ylim, '--', 'Color', [0.5667,0.5667,0.5667], 'LineWidth', 2)
ax.LineWidth = 2;

% save this figure
saveas(gcf, fullfile(filepath, [subID, '_PVT_plot.png']))
close all

end

%% Helper function
function [ data, fa, miss, temp_copy ] = readPVTtxt( filename )
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);
% Specify range and delimiter
opts.DataLines = [2, inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Trial", "RT", "ISI"];
opts.VariableTypes = ["string", "string", "string"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, ["Trial", "RT", "ISI"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Trial", "RT", "ISI"], "EmptyFieldRule", "auto");
% Import the data
temp = readmatrix(filename, opts);
% Clear temporary variables
clear opts
temp_copy = temp;

% Find the line reporting false alarms and misses
report_line = 0;
for ii = 1:size(temp, 1)
    if strcmpi(temp(ii,1), 'False Alarms')
        report_line = ii;
        assert(strcmpi(temp(ii,2), 'Misses'), 'Something is wrong with this .txt file.')
        fa = str2num(temp(ii+1,1)); 
        assert(~isempty(fa), 'No valid false alarm number extracted.')
        if ismissing(temp(ii+1,3))
            miss = str2num(temp(ii+1,2));
        else
            miss = str2num(temp(ii+1,3));
        end
        assert(~isempty(miss), 'No valid miss number extracted.')
    end
end
% get rid of the report lines 
temp(report_line:end, :) = [];

data = [];
data(:,1) = arrayfun(@(x) str2num(x), temp(:,1));
data(:,2) = arrayfun(@(x) str2num(x), temp(:,2));
data(:,3) = arrayfun(@(x) str2num(x), temp(:,3));

assert(all(~isnan(data), 'all'), 'data extracted is not valid!')

end



