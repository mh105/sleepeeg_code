function varargout = SleepEEG_loadedf(varargin)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - LOADEDF**
% 
% - this is a function used to load EDF/EDF+ file exported from the Natus
% clinical sleep EEG system at the MGH Sleep Center
% 
% This function is based on the blockEdfLoad() function from Dr. Dennis A.
% Dean, II, Ph.D. Modifications are made in order to accommodate the 'EDF
% Annotations' channel exported from the clinical system. The rest of the
% header, signal header, and signal vector specfications of the exported
% EDF/EDF+ file are the same as the original EDF format (Kemp et al. 1992).
%
% However, the content largely follows the 12 specifications recommended by
% EDF+, although some of the channel fields like "transducer_type" and
% "prefiltering" are not reported.
% 
% Check https://www.edfplus.info/specs/edfplus.html for full specs.
%
% ##### Function Prototypes:
%
%                                header = SleepEEG_loadedf(edfFN)
%                [header, signalHeader] = SleepEEG_loadedf(edfFN)
%    [header, signalHeader, signalCell] = SleepEEG_loadedf(edfFN)
%    [header, signalHeader, signalCell] = SleepEEG_loadedf(edfFN, signalLabels)
%    [header, signalHeader, signalCell] = SleepEEG_loadedf(edfFN, signalLabels, epochs)
%
% ##### Additional output (4th varargout):
%
%         --- annotation   this is a structure containing the onset time,
%                          duration, and content of various annotations
%                          made by sleep techs during overnight recording.
%
% ##### To load just the annotations from the EDF+ file:
%
%    [header, signalHeader, signalCell, annotation] = SleepEEG_loadedf(edfFN, {})
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% blockEdfLoad() Load EDF with memory block reads.
% Function inputs an EDF file text string and returns the header, signal
% header and each of the signals.
%
% The loader is designed to load the EDF file described in:
%
%    Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen and John Gade
%    "A simple format for exchange of digitized polygraphic recordings"
%    Electroencephalography and Clinical Neurophysiology, 82 (1992):
%    391-393.%
%
% An online description of the EDF format can be found at:
% http://www.edfplus.info/
%
% Our EDF tools can be found at:
%                  http://sleep.partners.org/edf/
%
% Updated versions will be posted:
%                  http://https://github.com/DennisDean/
%
% Additional Documentation and examples will be posted:
%                  http://sleepdata.org/tools
%
% Requirements:    Self contained, no external references
% MATLAB Version:  Requires R14 or newer, Tested with MATLAB 7.14.0.739
%
% Input (VARARGIN):
%           edfFN : File text string
%    signalLabels : Cell array of signal labels to return (optional)
%
% Output (VARARGOUT):
%          header : A structure containing variables for each header entry
%    signalHeader : A structured array containing signal information,
%                   for each structure present in the data
%      signalCell : A cell array that contains the data for each signal
%
% Output Structures:
%    header:
%       edf_ver
%       patient_id
%       local_rec_id
%       recording_startdate
%       recording_starttime
%       num_header_bytes
%       reserve_1
%       num_data_records
%       data_record_duration
%       num_signals
%    signalHeader (structured array with entry for each signal):
%       signal_labels
%       tranducer_type
%       physical_dimension
%       physical_min
%       physical_max
%       digital_min
%       digital_max
%       prefiltering
%       samples_in_record
%       reserve_2
%
% Examples:
%
%  Get EDF header information
%
%    edfFn3 = 'file.edf';
%    header = blockEdfLoad(edfFn3);
%
%    edfFn3 = 'file.edf';
%    [header signalHeader] = blockEdfLoad(edfFn3);
%
%
%  Load Signals
%
%    edfFn3 = 'file.edf';
%    [header signalHeader signalCell] = blockEdfLoad(edfFn3);
%
%    edfFn3 = 'file.edf';
%    signalLabels = {'Pleth', 'EKG-R-EKG-L', 'Abdominal Resp'};
%    [header signalHeader signalCell] = blockEdfLoad(edfFn3, signalLabels);
%
%    epochs = [1 2];  % Load first through second epoch
%    signalLabels = {'Pleth', 'Abdominal Resp', 'EKG-R-EKG-L'};
%    [header signalHeader signalCell] = ...
%         blockEdfLoad(edfFn3, signalLabels, epochs);
%
%
% Version: 0.1.21
%
% ---------------------------------------------
% Dennis A. Dean, II, Ph.D
%
% Program for Sleep and Cardiovascular Medicine
% Brigam and Women's Hospital
% Harvard Medical School
% 221 Longwood Ave
% Boston, MA  02149
%
% File created: October 23, 2012
% Last updated: January 23, 2014
%
% Copyright © [2012] The Brigham and Women's Hospital, Inc.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%

%%
%------------------------------------------------------------ Process input
% Opertion Flags
RETURN_PHYSICAL_VALUES = 1;

% Defaults for optional parameters
signalLabels = {};      % Labels of signals to return
epochs = [];            % Start and end epoch to return

% Process input
if nargin == 1
    edfFN = varargin{1};
    signalLabels = {};
elseif nargin == 2 && nargout >= 2
    edfFN = varargin{1};
    signalLabels = varargin{2};
elseif nargin == 3 && nargout >= 3
    edfFN = varargin{1};
    signalLabels = varargin{2};
    epochs = varargin{3};
else
    % Echo supported function prototypes to console
    fprintf('header = blockEdfLoad(edfFN)\n');
    fprintf('[header, signalHeader] = blockEdfLoad(edfFN)\n');
    fprintf('[header, signalHeader, signalCell] = blockEdfLoad(edfFN)\n');
    fprintf('[header, signalHeader, signalCell] = blockEdfLoad(edfFN, signalLabels)\n');
    fprintf('[header, signalHeader, signalCell] = blockEdfLoad(edfFN, signalLabels, epochs)\n');
    
    % Call MATLAB error function
    error('Function prototype not valid');
end

%-------------------------------------------------------------- Input check
% Check that first argument is a string
if ~ischar(edfFN)
    msg = ('First argument is not a string.');
    error(msg);
end
% Check that second argument is a cell string
if ~iscellstr(signalLabels)
    msg = ('Second argument is not a valid text string.');
    error(msg);
end
% Check that third argument is valid with two Epoch times
if and(nargin ==3, length(epochs)~=2)
    msg = ('Specify epochs = [Start_Epoch End_Epoch]');
    error(msg);
end

%---------------------------------------------------  Load File Information
% Open file for reading
[fid, msg] = fopen(edfFN);

% Proceed if file is valid
if fid <0
    % file id is not valid
    error(msg);
end

%%
%-------------------------------------------------------------- Load Header
try
    % Load header information in one call
    edfHeaderSize = 256;
    [A, count] = fread(fid, edfHeaderSize);
    assert(count == edfHeaderSize, ['Not all ', num2str(edfHeaderSize) ' are successfully loaded!'])
catch exception
    msg = 'File load error. Cannot read header information. Check available memory.';
    error(msg);
end

%----------------------------------------------------- Process Header Block
% These specifications are hard-coded by the EDF format.
% See https://www.edfplus.info/specs/edf.html for descriptions.

% Create array/cells to create struct with loop
headerVariables = {...
    'edf_ver';            'patient_id';         'local_rec_id'; ...
    'recording_startdate';'recording_starttime';'num_header_bytes'; ...
    'reserve_1';          'num_data_records';   'data_record_duration';...
    'num_signals'};
headerVariablesConF = {...
    @strtrim;   @strtrim;   @strtrim; ...
    @strtrim;   @strtrim;   @str2num; ...
    @strtrim;   @str2num;   @str2num;...
    @str2num};
headerVariableSize = [8; 80; 80; 8; 8; 8; 44; 8; 8; 4];
headerVarLoc = vertcat([0],cumsum(headerVariableSize));
headerSize = sum(headerVariableSize);

% Create Header Structure
header = struct();
for h = 1:length(headerVariables)
    conF = headerVariablesConF{h};
    value = conF(char(A(headerVarLoc(h)+1:headerVarLoc(h+1))'));
    header.(headerVariables{h}) = value;
end

% End Header Load section

%%
%------------------------------------------------------- Load Signal Header
if nargout >= 2
    try
        % Load signal header into memory in one load
        edfSignalHeaderSize = header.num_header_bytes - headerSize;
        [A, count] = fread(fid, edfSignalHeaderSize);
        assert(count == edfSignalHeaderSize, ['Not all ', num2str(edfSignalHeaderSize) ' are successfully loaded!'])
    catch exception
        msg = 'File load error. Cannot read signal header information. Check available memory.';
        error(msg);
    end
    
    %------------------------------------------ Process Signal Header Block
    % These specifications are hard-coded by the EDF format.
    % See https://www.edfplus.info/specs/edf.html for descriptions.
    
    % Create array/cells to create struct with loop
    signalHeaderVar = {...
        'signal_labels'; 'tranducer_type'; 'physical_dimension'; ...
        'physical_min'; 'physical_max'; 'digital_min'; ...
        'digital_max'; 'prefiltering'; 'samples_in_record'; ...
        'reserve_2' };
    signalHeaderVarConvF = {...
        @strtrim; @strtrim; @strtrim; ...
        @str2num; @str2num; @str2num; ...
        @str2num; @strtrim; @str2num; ...
        @strtrim };
    num_signal_header_vars = length(signalHeaderVar);
    num_signals = header.num_signals;
    signalHeaderVarSize = [16; 80; 8; 8; 8; 8; 8; 80; 8; 32];
    signalHeaderBlockSize = sum(signalHeaderVarSize)*num_signals;
    signalHeaderVarLoc = vertcat([0],cumsum(signalHeaderVarSize*num_signals));
    signalHeaderRecordSize = sum(signalHeaderVarSize);
    
    % Create Signal Header Struct
    signalHeader = struct(...
        'signal_labels', {},'tranducer_type', {},'physical_dimension', {}, ...
        'physical_min', {},'physical_max', {},'digital_min', {},...
        'digital_max', {},'prefiltering', {},'samples_in_record', {},...
        'reserve_2', {});
    
    % Get each signal header variable
    for v = 1:num_signal_header_vars
        varBlock = A(signalHeaderVarLoc(v)+1:signalHeaderVarLoc(v+1))';
        varSize = signalHeaderVarSize(v);
        conF = signalHeaderVarConvF{v};
        for s = 1:num_signals
            varStart = varSize*(s-1)+1;
            varEnd = varSize*s;
            value = conF(char(varBlock(varStart:varEnd)));
            signalHeader(s).(signalHeaderVar{v}) = value;
        end
    end
    
end % End Signal Header Load Section

%%
%-------------------------------------------------------- Load Signal Block
if nargin >= 2 || nargout >=3
    % Read digital values to the end of the file
    
    % Added functionality to be compatible with EDF+ feature in the
    % exported EDF/EDF+ file from Natus clinical sleep EEG system.
    % Specifically, we look for a signal channel with the label 'EDF
    % Annotations' so that we can load it differently from int16 encoding.
    annot = arrayfun(...
        @(x)strcmp(signalHeader(x).signal_labels,'EDF Annotations'),...
        1:header.num_signals, 'UniformOutput', false);
    annotidx = find(cell2mat(annot)==1);
    
    try
        % Set default error message
        errMsg = 'File load error. Cannot read signal values. Check available memory.';
        
        % Load strategy is dependent on input
        if nargin == 1
            % Load entire file
            [A, ~] = fread(fid, 'int16');
            
            edfSignalSizes = arrayfun(...
                @(x)signalHeader(x).samples_in_record, [1:header.num_signals]);
            edfRecordSize = sum(edfSignalSizes);
        else
            % Get signal label information
            edfSignalLabels = arrayfun(...
                @(x)signalHeader(x).signal_labels, [1:header.num_signals],...
                'UniformOutput', false);
            signalIndexes = arrayfun(...
                @(x)find(strcmp(x,edfSignalLabels)), signalLabels,...
                'UniformOutput', false);
            
            % Check that specified signals are present
            signalIndexesCheck = cellfun(...
                @(x)~isempty(x), signalIndexes, 'UniformOutput', false);
            signalIndexesCheck = int16(cell2mat(signalIndexesCheck));
            if sum(signalIndexesCheck) == length(signalIndexes)
                % Indices are specified
                signalIndexes = cell2mat(signalIndexes);
            else
                % Couldn't find at least one signal label
                errMsg = 'Could not identify signal label. Please check!';
                error(errMsg);
            end
            
            edfSignalSizes = arrayfun(...
                @(x)signalHeader(x).samples_in_record, [1:header.num_signals]);
            edfRecordSize = sum(edfSignalSizes);
            
            % Identify memory locations to record
            endLocs = cumsum(edfSignalSizes)';
            startLocs = [1;endLocs(1:end-1)+1];
            signalLocs = [];
            for s = signalIndexes
                signalLocs = [signalLocs; [startLocs(s):1:endLocs(s)]'];
            end
            sizeSignalLocs = length(signalLocs);
            
            % Load only required signals reduce memory calls
            loadedSignalMemory = header.num_data_records*...
                sum(edfSignalSizes(signalIndexes));
            A = zeros(loadedSignalMemory,1);
            for r = 1:header.num_data_records
                [a, count] = fread(fid, edfRecordSize, 'int16');
                assert(count == edfRecordSize, ['Not all ', num2str(edfRecordSize) ' are successfully loaded!'])
                A([1+sizeSignalLocs*(r-1):sizeSignalLocs*r]) = a(signalLocs);
            end
            
            % Reset global variables, which enable reshape functions to
            % work correctly
            header.num_signals = length(signalLabels);
            signalHeader = signalHeader(signalIndexes);
            num_signals = length(signalIndexes);
        end
        
    catch exception
        error(errMsg);
    end
    
    %------------------------------------------------- Process Signal Block
    % Get values to reshape block
    num_data_records = header.num_data_records;
    getSignalSamplesF = @(x)signalHeader(x).samples_in_record;
    signalSamplesPerRecord = arrayfun(getSignalSamplesF,[1:num_signals]);
    recordWidth = sum(signalSamplesPerRecord);
    
    % Reshape - Each row is a data record
    A = reshape(A, recordWidth, num_data_records)';
    
    % Create raw signal cell array
    signalCell = cell(1,num_signals);
    signalLocPerRow = horzcat([0],cumsum(signalSamplesPerRecord));
    for s = 1:num_signals
        % Get signal location
        signalRowWidth = signalSamplesPerRecord(s);
        signalRowStart = signalLocPerRow(s)+1;
        signaRowEnd = signalLocPerRow(s+1);
        
        % Create Signal
        signal = reshape(A(:,signalRowStart:signaRowEnd)',...
            signalRowWidth*num_data_records, 1);
        
        % Get scaling factors
        dig_min = signalHeader(s).digital_min;
        dig_max = signalHeader(s).digital_max;
        phy_min = signalHeader(s).physical_min;
        phy_max = signalHeader(s).physical_max;
        
        % Assign signal value
        value = signal;
        
        % Convert to physical units
        if RETURN_PHYSICAL_VALUES == 1
            % Convert from digital to physical values
            value = (signal-dig_min)/(dig_max-dig_min);
            value = value.*double(phy_max-phy_min)+phy_min;
        else
            fprintf('Digital to Physical conversion is NOT performned: %s\n',...
                edfFN);
        end
        
        signalCell{s} = value;
    end
    
end % End Signal Load Section

%%
%-------------------------------------------------- Extract Annotations
if nargout ==4
    % Rewind the fid pointer to the start of the file
    frewind(fid);
    % Burn out the header lines
    [~, ~] = fread(fid, header.num_header_bytes);
    % Now let's load the annotations into memory
    annotation = {};
    edfRecordSize = sum(edfSignalSizes);
    for r = 1:header.num_data_records
        % first burn everything before the annotation channel
        [~, ~] = fread(fid, sum(edfSignalSizes(1:annotidx-1)), 'int16');
        % then we load the annotation channel:
        % for some reason, Natus system encodes these annotations as
        % samples_in_record (ns) * 2-byte integers, but the precision of
        % the characters are actually encoded in uint8 (unsigned integer
        % in 8 bits, single-byte). As the result, there are actually 2
        % times ns bytes. We need to take this into account when using
        % fread.
        [a, count] = fread(fid, edfSignalSizes(annotidx)*2);
        assert(count == edfSignalSizes(annotidx)*2, ['Not all ', num2str(edfSignalSizes(annotidx)*2) ' are successfully loaded!'])
        % this should have exhausted one data record
        assert(edfRecordSize - sum(edfSignalSizes(1:annotidx)) == 0, 'Non-empty remaining samples in a data record.')
        % now let's parse the annotation strings into Time-stamped
        % Annotations Lists (TALs)
        tempstr = char(a');
        if contains(tempstr, char([20, 0]))
            tempstr = strsplit(tempstr, char([20, 0]));
            annotation(size(annotation,1)+1:size(annotation,1)+length(tempstr), 1) = tempstr';
        else
            error('Current annotation entry in the data record does not have char([20, 0])')
        end
    end
    
    % Now we need to organize the annotation cell into a human readable
    % form as a structure
    annotation = SleepEEG_annot(header, annotation);

end % End Annotation Load Section

%%
%------------------------------------------------------ Create return value
if nargout < 2
    varargout{1} = header;
elseif nargout == 2
    varargout{1} = header;
    varargout{2} = signalHeader;
elseif nargout >= 3
    % Check if a reduce signal set is requested
    if ~isempty(epochs)
        % Determine signal sampling rate
        signalSamples = arrayfun(...
            @(x)signalHeader(x).samples_in_record, [1:num_signals]);
        signalIndex = ones(num_signals, 1)*[epochs(1)-1 epochs(2)]*30;
        samplesPerSecond = (signalSamples/header.data_record_duration)';
        signalIndex = signalIndex .* [samplesPerSecond samplesPerSecond];
        signalIndex(:,1) = signalIndex(:,1)+1;
        
        % Redefine signals to include specified epochs
        signalIndex = int64(signalIndex);
        for s = 1:num_signals
            signal = signalCell{s};
            index = [signalIndex(s,1):signalIndex(s,2)];
            signalCell{s} = signal(index);
        end
    end
    
    % Create Output Structure
    varargout{1} = header;
    varargout{2} = signalHeader;
    varargout{3} = signalCell;
    
    if nargout ==4 % 4th output argument to return the annotation list
        varargout{4} = annotation;
    end 
    
end % End Return Value Function

%%
% Close file explicitly
if fid > 0
    fclose(fid);
end

end % End of SleepEEG_loadedf() function
