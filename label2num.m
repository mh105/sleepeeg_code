function [ channelnum ] = label2num(stringc, chanlocs)
    for ii = 1:length(chanlocs)
        if strcmp(chanlocs(ii).labels, stringc)
            channelnum = ii;
        end
    end
    assert(exist('channelnum', 'var')==1, 'No channel number found!')
end