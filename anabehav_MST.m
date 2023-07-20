function [ subData ] = anabehav_MST(subID, filepath)
%
% **ADSLEEPEEG_COGNITIVE TASK BEHAVIORAL ANALYSES FUNCTION - MST**
%
% - Motor Sequence Task (MST)
%
% - loads a set of 4 text files with MST performance saved by NBS
% Presentation experiments, plots basic behavioral performance, and saves a
% single .mat file with organized behavioral performance data
%
% **Need to add _lefty related codes! 
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%                           this function will automatically search for the
%                           4 MST task files in the filepath directory:
%                               - subID_night1_MST_V1_Test1.txt
%                               - subID_morning1_MST_V1_Test2.txt
%                               - subID_night2_MST_V2_Test1.txt
%                               - subID_morning2_MST_V2_Test2.txt
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
try 
    v1t1 = fullfile(filepath, [subID, '_night1_MST_V1_Test1.txt']);
    v1t2 = fullfile(filepath, [subID, '_morning1_MST_V1_Test2.txt']);
    v2t1 = fullfile(filepath, [subID, '_night2_MST_V2_Test1.txt']);
    v2t2 = fullfile(filepath, [subID, '_morning2_MST_V2_Test2.txt']);
    assert(isfile(v1t1))
    assert(isfile(v1t2))
    assert(isfile(v2t1))
    assert(isfile(v2t2))
catch 
    v1t1 = fullfile(filepath, [subID, '_night2_MST_V1_Test1.txt']);
    v1t2 = fullfile(filepath, [subID, '_morning2_MST_V1_Test2.txt']);
    v2t1 = fullfile(filepath, [subID, '_night1_MST_V2_Test1.txt']);
    v2t2 = fullfile(filepath, [subID, '_morning1_MST_V2_Test2.txt']);
end

assert(isfile(v1t1), [subID, '_MST_V1_Test1.txt file is missing from the directory!'])
assert(isfile(v1t2), [subID, '_MST_V1_Test2.txt file is missing from the directory!'])
assert(isfile(v2t1), [subID, '_MST_V2_Test1.txt file is missing from the directory!'])
assert(isfile(v2t2), [subID, '_MST_V2_Test2.txt file is missing from the directory!'])

%% Load data files
% Start loading...

sID = []; % Subject ID
ExpVersion = []; % Two versions: 1 and 2
SessionNum = []; % Two test sessions: 1 and 2
NightNum = []; % Night number, 1 = night1/morning1, 2 = night2/morning2
blockNum = []; % Block number: 1-12
seqNum = []; % Sequence number - the order of completed sequence in a block
meanRT = []; % mean Response time in msec for completing each sequence
rawtxt = {}; % store of raw data imported from the 4 text files

for v = 1:2 % Version 1 and 2
    if v==1
        nightn = str2num(v1t1(strfind(v1t1, 'night')+5));
    else
        nightn = str2num(v2t1(strfind(v2t1, 'night')+5));
    end
    for s = 1:2 % Test Session 1 and 2
        if v==1 && s==1
            textfn = v1t1;
        elseif v==1 && s==2
            textfn = v1t2;
        elseif v==2 && s==1
            textfn = v2t1;
        elseif v==2 && s==2
            textfn = v2t2;
        end
        
        % specify correct sequences depending on the version
        if v==1; corrSeq = [4 1 3 2 4];else; corrSeq = [2 4 2 3 1];end
        
        % load from text file
        delimiterIn = '\t';
        headerlinesIn = 1;
        A = importdata(textfn, delimiterIn, headerlinesIn);
        data = A.data;
        A.textfn = textfn;
        A.version = v;
        A.testsession = s;
        rawtxt{length(rawtxt)+1} = A;
        
        % Find correct sequences completed within each of 12 blocks
        for i = 1:length(unique(data(:,1)))
            temp_data = data(data(:,1)==i, :); % look within a block
            Index  = strfind(temp_data(:,3)', corrSeq);
            temp_data(2:end,4) = diff(temp_data(:,4));
            
            tempsID = string;
            tempsID(1:length(Index),:) = subID;
            sID = [sID; tempsID];
            ExpVersion = [ExpVersion; ones(length(Index), 1).*v];
            SessionNum = [SessionNum; ones(length(Index), 1).*s];
            
            NightNum = [NightNum; ones(length(Index), 1).*nightn];
            blockNum = [blockNum; ones(length(Index), 1).*i];
            seqNum = [seqNum; [1:length(Index)]']; %#ok<NBRAK>
            
            for k = 1:length(Index) % Iterate through the sequences completed
                % compute mean RT for each completed sequence with 5 presses
                meanRT = [meanRT; mean(temp_data(Index(k):Index(k)+4, 4))]; 
            end   
        end
    end
end

% Organizing into a structure variable
subData = struct;
subData.subID = sID;
subData.ExpVersion = ExpVersion;
subData.SessionNum = SessionNum;
subData.NightNum = NightNum;
subData.blockNum = blockNum;
subData.seqNum = seqNum;
subData.meanRT = meanRT;
subData.rawtxt = rawtxt;

% save the structure 
save(fullfile(filepath, [subID, '_MST_data']), 'subData')

%% Behavioral analysis

% RT measures
seqCount = [];
meanRT = [];

% ExpVersion = []; % Two versions: 1 and 2
% SessionNum = []; % Two test sessions: 1 and 2
% blockNum = []; % Block number: 1-12
for m = 1:length(unique(subData.ExpVersion))
    for j = 1:length(unique(subData.SessionNum))
        for k = 1:length(unique(subData.blockNum))
            sel = subData.ExpVersion == m & subData.SessionNum == j & subData.blockNum == k;
            seqCount(m,j,k) = length(subData.seqNum(sel));
            % averaging across sequences within a block
            meanRT(m,j,k) = mean(subData.meanRT(sel)); 
        end
    end
end

%Plot overall Performance averaged across all sessions
figure
set(gcf,'units', 'pixels','Position',[0 0 1000 800])
% averaging across the two versions

subplot(2,1,1)
hold on
presleep_count = squeeze(mean(seqCount(:,1,:), 1));
postsleep_count = squeeze(mean(seqCount(:,2,:), 1));
s1 = scatter(1:length(presleep_count), presleep_count, 100, 'filled');
s2 = scatter([1:length(postsleep_count)]+length(presleep_count)+1, postsleep_count, 100, 'filled');
title('MST: Completed Sequences in 30-sec Blocks')
ylabel('Sequence Count')
xticks([6, 13, 20])
xticklabels({'Pre-sleep Block 1-12', 'SLEEP', 'Post-sleep Block 1-12'})
set(gca, 'FontSize', 20)
legend([s1, s2], {'Pre-sleep', 'Post-sleep'}, 'Location', 'best')
drawnow
ylimit = ylim;
rectangle('Position',[length(presleep_count)+0.5,ylimit(1),1,ylimit(2)-ylimit(1)],'FaceColor',[0.6 0.6 0.9],'EdgeColor','none',...
    'LineWidth',0.5);

subplot(2,1,2)
hold on
presleep_count = squeeze(mean(meanRT(:,1,:), 1));
postsleep_count = squeeze(mean(meanRT(:,2,:), 1));
s1 = scatter(1:length(presleep_count), presleep_count, 100, 'filled');
s2 = scatter([1:length(postsleep_count)]+length(presleep_count)+1, postsleep_count, 100, 'filled');
title('MST: Mean RT of Completed Sequences in 30-sec Blocks')
ylabel('mean RT (ms)')
xticks([6, 13, 20])
xticklabels({'Pre-sleep Block 1-12', 'SLEEP', 'Post-sleep Block 1-12'})
set(gca, 'FontSize', 20)
legend([s1, s2], {'Pre-sleep', 'Post-sleep'}, 'Location', 'best')
drawnow 
ylimit = ylim;
rectangle('Position',[length(presleep_count)+0.5,ylimit(1),1,ylimit(2)-ylimit(1)],'FaceColor',[0.6 0.6 0.9],'EdgeColor','none',...
    'LineWidth',0.5);

% save this figure
saveas(gcf, fullfile(filepath, [subID, '_MST_plot.png']))
close all

end

