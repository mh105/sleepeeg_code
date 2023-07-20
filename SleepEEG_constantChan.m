function [ Constantidx, offsetlist, EEG ] = SleepEEG_constantChan(EEG, cthresh, setzero, ChanNum)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - CONSTANTCHAN**
%
% - Used to zero-center channels recording mostly constant values
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - EEG:          data structure containing the EEG data and
%                           impedance values.
%
%           - cthresh:      percentage threshold for identifying as a
%                           constant channel.
%                           default: 0.95
%
%           - setzero:      whether to set the constant channels to zero.
%                           default: false
%
%           - ChanNum:      channel numbers to check for constant values.
%                           default: [1]
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - Constantidx:  a vector of double type containing the channel
%                           number of electrodes detected to have constant
%                           values throughout the recording.
%
%           - offsetlist:   a vector of double type with the same size as
%                           the first input Constantidx containing the
%                           offset values of the constant channels from
%                           zero.
%
%           - EEG:          an updated EEG structure if setzero was true
%                           then EEG.data for the constant channels are set
%                           to be be a constant zero recording.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('---------------------------------')
disp([mHead 'SleepEEG_constantChan()']);
disp('---------------------------------')

%% Set default threshold
if nargin < 2
    cthresh = 0.95;
    setzero = false;
elseif nargin < 3
    setzero = false;
end

if ~exist('ChanNum', 'var')
    % since we often don't gel LM electrode, we check just this channel.
    ChanNum = 1;
end

%% Constant electrode detection
signallength = EEG.pnts;
threshlength = floor(signallength * cthresh);

Constantidx = [];
offsetlist = [];

for ii = ChanNum
    disp([mHead 'Checking for constant values in channel: ', num2str(ii)])
    dxChan = [0, diff(EEG.data(ii,:))];
    n0length = sum(dxChan == 0);
    if n0length >= threshlength
        % this electrode has mostly constant values
        Constantidx = [Constantidx, ii]; %#ok<*AGROW>
        
        % find the offset from zero
        offset = mode(EEG.data(ii,:));
        offsetlist = [offsetlist, offset];
        
        if setzero
            % set this channel to a constant zero for all timepoints
            disp([mHead 'Setting channel ', num2str(ii), ' to zero.'])
            EEG.data(ii,:) = zeros(size(EEG.data(ii,:)));
        end
    end
end
disp([mHead 'Finished checking for constant values in channels.'])

if isempty(Constantidx)
    disp([mHead 'No constant channel was found.'])
elseif setzero && ~isempty(Constantidx)
    disp([mHead 'Constant channels are set to zero for all timepoints: ', num2str(Constantidx)]);
end

end
