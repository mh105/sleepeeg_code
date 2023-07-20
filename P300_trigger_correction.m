% Correct the button response trigger code in P300 to disambiguate it from
% the oddball trial trigger code 
%
% Last edit: Alex 02/24/2022

% First find the block start and end latencies

block1_start = 0;
block1_end = 0;
block2_start = 0;
block2_end = 0;

for ii = 1:length(EEG.event)
    if strcmp(EEG.event(ii).type, '22')
        if block1_start == 0
            block1_start = EEG.event(ii).latency;
        else
            block2_start = EEG.event(ii).latency;
        end
    elseif strcmp(EEG.event(ii).type, '23')
        if block1_end == 0
            block1_end = EEG.event(ii).latency;
        else
            block2_end = EEG.event(ii).latency;
        end
    end
end

% Then find all the '3' trigger codes during testing blocks 
trigger_3_blocks = [];

for ii = 1:length(EEG.event)
    if strcmp(EEG.event(ii).type, '3')
        latency = EEG.event(ii).latency;
        block_index = size(trigger_3_blocks,1) + 1;
        if latency > block1_start && latency < block1_end
            trigger_3_blocks(block_index, 1) = ii;
        elseif latency > block2_start && latency < block2_end
            trigger_3_blocks(block_index, 1) = ii;
        else
            EEG.event(ii).type = '11';
        end
    end
end

assert(length(trigger_3_blocks) == 80, 'Not 100% accuracy, you need to check manually.')

for l = 2:2:size(trigger_3_blocks, 1)
    EEG.event(trigger_3_blocks(l)).type = '11';
end
