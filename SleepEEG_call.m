function [ EEG1 ] = SleepEEG_call(subID, fnsuffix, varargin)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - CALL**
%&#x1F536;
%
% - this is the master function to be called in the terminal window. One can
% call this function with subject ID and file suffix with various processing
% arguments that will be applied in an appropriate order
%
% Usage: [ EEG1 ] = SleepEEG_call(subID, fnsuffix, '<flag#1>',<arg#1>...'<flag#n>',<arg#n>);
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Declared Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           EEG .cnt file, usually "night1(2)_Sleep". if specified as
%                           [] empty, then all .cnt files under the subID
%                           directory will be processed.
%
%                           ***the final file name is in the form: subID_fnsuffix.cnt***
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Optional Inputs:
%
%           - []:           When SleepEEG_call is called without varargin
%                           nor fnsuffix, only 'report' is true, all other
%                           commands are set to false.
%
%           - 'report':     To report ALL DATA files collected for the
%                           subject, and prints in command window the types
%                           of preprocessing that have been completed for
%                           the subject's .cnt files.
%                           default: false
%
%           - 'verbose':    whether print messages during processing.
%                           default: true
%
%           - 'downsample': downsample to the a specified rate in Hz.
%                           default: [true, 500 (Hz)]
%
%           - 'dead':       whether to detect dead electrodes.
%                           default: true
%
%           - 'interp':     whether to interpolate dead electrodes.
%                           default: false
%
%           - 'reref':      whether re-reference to common average reference.
%                           options are:
%                               - AR (common average)
%                               - Z3 (recording reference at Z3)
%                               - [CH] (arbitrary channel number as reference)
%                               - left mastoid
%                               - linked mastoid
%                               - REST (reference electrode standardization
%                                       technique)
%                               - contral mastoid
%                               - LP (Laplacian reference based on duke layout)
%                           default: {false, 'AR'}
%
%           - 'oversave':   whether to overwrite existing .set files and
%                           save the EEG structure again onto disk.
%                           default: false
%
%           - 'impedance':  whether check impedance values at the beginning
%                           and the end of EEG recording.
%                           default: true
%
%           - 'bridge':     whether to do bridging analysis.
%                           default: true
%
%           - 'Pick2Spec':  channel number to plot the time traces,
%                           locations, and spectrograms of 1 anterior and 1
%                           posterior electrodes.
%                           default: [25, 84]
%
%           - 'originFs':   channel number to do spectrogram for two
%                           electrodes in original sampling rate without
%                           downsampling.
%                           default: []
%
%           - project:      an optional string to specify project name for
%                           path configuration.
%                           default: ''
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

% Set warning of graphical display errors in command window to off
warning('off', 'MATLAB:callback:error')

%% Configure input arguments
% Construct commandDict structure:
commandHeadings = {'report','verbose','downsample','dead','interp','reref','oversave',...
    'impedance','bridge','Pick2Spec','originFs','project'}; % command field names

% Some pre-set command setting keywords:
if any(strcmpi('all', varargin))
    disp([mHead, 'Preset command "all" used, additional arguments will be updated.']);
    commandDefaults = {true,true,[true,500],true,false,{false,'AR'},false,true,...
        true,[25,84],[25,84],''};
elseif any(strcmpi('checkgel', varargin))
    disp([mHead, 'Preset command "checkgel" used, additional arguments will be updated.']);
    commandDefaults = {false,true,[true,500],true,false,{false,'AR'},false,true,...
        true,[25,84],[],''};
elseif any(strcmpi('none', varargin))
    disp([mHead, 'Preset command "none" used, additional arguments will be updated.']);
    commandDefaults = {false,true,[false,500],false,false,{false,'AR'},false,false,...
        false,[],[],''};
elseif any(strcmpi('spectro', varargin))
    disp([mHead, 'Preset command "spectro" used, additional arguments will be updated.']);
    commandDefaults = {false,true,[true,500],false,false,{false,'AR'},false,false,...
        false,[25,84],[],''};
elseif nargin < 2 % if only subID provided -> call SleepEEG_report.m
    disp([mHead, 'Only subID provided, will call SleepEEG_report().']);
    commandDefaults = {true,false,[false,500],false,false,{false,'AR'},false,false,...
        false,[],[],''}; fnsuffix = [];
else % default settings
    disp([mHead, 'Pulling default commands, additional arguments will be updated.']);
    commandDefaults = {false,true,[true,500],true,false,{false,'AR'},false,true,...
        true,[25,84],[],''};
end

% Update commandDict values according to varargin
for i = 1:length(commandHeadings)
    FlagIndex = find(strcmpi(commandHeadings{i}, varargin)==1);
    assert(length(FlagIndex) <= 1,'Only one %s value can be entered as an input.', commandHeadings{i})
    if ~isempty(FlagIndex)
        commandDefaults{i} = varargin{FlagIndex+1};
    end
end

% Create a dictionary structure storing the command argument variables
commandDict = cell2struct(commandDefaults, commandHeadings, 2);

clear commandHeadings commandDefaults FlagIndex

%% Recursion decision point: base case or loop through all available files?
if isempty(fnsuffix) && nargin > 2 % if fnsuffix is not specified, we process all available .cnt files for the subject
    [~, ~, fnsuffixlist, ~] = SleepEEG_configDir(subID, [], commandDict.verbose, commandDict.project);

    tstart = tic;
    for i = 1:length(fnsuffixlist)
        SleepEEG_call(subID, fnsuffixlist{i}, varargin{:});
    end
    if commandDict.verbose
        disp([mHead, 'Total time taken to analyze subject <<<', subID ,'>>>...']);
        fprintf(mSpace)
        toc(tstart)
    end

else % otherwise use fnsuffix to select one .cnt file

    %% Define directories of codes and data folders
    [dataDir, datafn, fileID, outputDir] = SleepEEG_configDir(subID, fnsuffix, commandDict.verbose, commandDict.project);

    %% Log command window outputs
    logDir = fullfile(dataDir, subID, 'log');
    if ~exist(logDir, 'dir')
        mkdir(logDir) % make the directory if non-existant
    end

    startdatetime = datetime(datetime);
    startdatetime = strrep(startdatetime,':','_');

    if ~isempty(datafn)
        diary(fullfile(dataDir, subID, 'log', [fileID, '_log_', startdatetime, '.txt']))
    else
        diary(fullfile(dataDir, subID, 'log', ['_log_', startdatetime, '.txt']))
    end

    disp([mHead, 'Destination folder for log file:'])
    fprintf(mSpace)
    disp(logDir)
    disp(' ')

    %% Display the list of commands to be executed
    disp('-------------------------')
    disp([mHead 'SleepEEG_call()']);
    disp('-------------------------')
    disp([mHead 'Call commands to be executed:']);
    disp(commandDict)
    if commandDict.verbose; procstart = tic; end

    %% ----%%----%%----%%----%%----%%---
    %  ----%%----%%----%%----%%----%%---
    %
    %
    %     Start of Command Execution
    %
    %
    %  ----%%----%%----%%----%%----%%---
    %  ----%%----%%----%%----%%----%%---

    %% "report"
    %           - 'report':     To report ALL DATA files collected for the
    %                           subject, and prints in command window the types
    %                           of preprocessing that have been completed for
    %                           the subject's .cnt files.
    %                           default: true
    if commandDict.report
        SleepEEG_report(subID, commandDict.project)
    end

    %% "verbose"
    %           - 'verbose':    whether print messages during processing.
    %                           default: true
    if commandDict.verbose
        warning('on','all')
        disp([mHead 'Verbose is set to true.']);
    else
        warning('off','all')
        disp([mHead 'Verbose is set to false.']);
    end
    fprintf(mHead); warning

    %% "downsample"
    %           - 'downsample': downsample to the a specified rate in Hz.
    %                           default: [true, 500 (Hz)]
    if ~isempty(datafn)
        %   When downsample flag is set to false, SleepEEG_downsample will
        %   just load in the .cnt file and return it without downsampling.
        EEG1 = SleepEEG_downsample(subID, fnsuffix, commandDict.downsample,...
            commandDict.verbose, commandDict.oversave, commandDict.project,...
            dataDir, datafn);
    end

    %% 'dead'
    %           - 'dead':       whether to detect dead electrodes.
    %                           default: true
    if commandDict.dead
        % The recording reference electrode is removed from the Deadidx
        Deadidx = SleepEEG_dead(subID, fnsuffix, commandDict.project, EEG1, fileID, outputDir);
    end

    %% 'interp'
    %           - 'interp':     whether to interpolate dead electrodes.
    %                           default: false
    if commandDict.interp
        % Before we do re-referencing, we need to interpolate and recover
        % the dead electrodes except the original reference electrode
        disp([mHead, 'Dead electrodes interpolated using griddata function with method [v4].'])
        EEG1 = pop_interp(EEG1, Deadidx, 'v4');
    end

    %% "reref"
    %           - 'reref':      whether re-reference to common average reference.
    %                           options are:
    %                               - AR (common average)
    %                               - Z3 (recording reference at Z3)
    %                               - [CH] (arbitrary channel number as reference)
    %                               - left mastoid
    %                               - linked mastoid
    %                               - REST (reference electrode standardization
    %                                       technique)
    %                               - contral mastoid
    %                               - LP (Laplacian reference based on duke layout)
    %                           default: {false, 'AR'}
    if commandDict.reref{1}
        % channel 129 is EOG, excluded from common average so it is kept
        % as unipolar referenced to the recording reference electrode (Z3)
        EEG1 = ANT_interface_reref([], commandDict.reref{2}, [], 129, EEG1, commandDict.verbose);
    else
        disp([mHead, 'No re-referencing was applied.'])
    end

    %% 'oversave'
    %           - 'oversave':   whether to overwrite existing .set files and
    %                           save the EEG structure again onto disk. Use
    %                           this flag to save re-referenced or
    %                           interpolated data by setting it to true.
    %                           default: false
    if commandDict.oversave
        % Specify .set file name
        savefn = [subID, '_', fnsuffix];
        if commandDict.downsample(1)
            savefn = [savefn, '_ds', num2str(commandDict.downsample(2))];
        end
        savefn = [savefn, '_', EEG1.refscheme];
        % Specify .set file folder
        subfolder = 'set';
        filepath = fullfile(dataDir, subID, subfolder);
        % only save if oversave is true
        EEG1 = ANT_interface_saveset(EEG1, savefn, filepath, commandDict.verbose);
    end

    %% ----%%----%%----%%----%%----%%---%%----%%---%%----%%---%%----%%---
    % All following commands are analysis commands that will be applied
    % on EEG structure when available. However, no change will be made to
    % the EEG structure, so no saving will be attempted after this point.
    %  ----%%----%%----%%----%%----%%---%%----%%---%%----%%---%%----%%---

    % Report output directory
    if commandDict.verbose
        disp(' ')
        disp([mHead, 'Analyzing the loaded data...'])
        disp(' ')

        disp([mHead, 'Data analysis output directory:'])
        fprintf(mSpace);
        disp(outputDir)
        disp(' ')
    end

    %% 'impedance'
    %           - 'impedance':  whether check impedance values at the beginning
    %                           and the end of EEG recording.
    %                           default: true
    if commandDict.impedance
        SleepEEG_impcheck(subID, fnsuffix, commandDict.project, EEG1, fileID, outputDir);
    end

    %% 'bridge'
    %           - 'bridge':     whether to do bridging analysis.
    %                           default: true
    if commandDict.bridge
        SleepEEG_bridge(subID, fnsuffix, commandDict.project, EEG1, fileID, outputDir);
    end

    %% 'Pick2Spec'
    %           - 'Pick2Spec':  channel number to plot the time traces,
    %                           locations, and spectrograms of 1 anterior and 1
    %                           posterior electrodes.
    %                           default: [25, 84]
    if ~isempty(commandDict.Pick2Spec)
        SleepEEG_AntPosSpec(subID, fnsuffix, commandDict.project, EEG1, commandDict.Pick2Spec, fileID, outputDir);
    end

    %% 'originFs'
    %           - 'originFs':   channel number to do spectrogram for two
    %                           electrodes in original sampling rate without
    %                           downsampling.
    %                           default: []
    if ~isempty(commandDict.originFs)

        % This is a debugging tool to confirm that no artifact is introduced
        % in the downsampling process.

        % Extract two channels
        [~, srate] = SleepEEG_extractChan(subID, fnsuffix, commandDict.originFs, commandDict.project, dataDir, datafn, fileID);

        % Reconstruct the two channels
        channelfn = SleepEEG_reconChan(subID, fnsuffix, commandDict.originFs, commandDict.project, dataDir, datafn, fileID);

        % Load in the channel data
        data = [];
        for i = 1:length(commandDict.originFs)
            load(channelfn{i}, 'channel_data')
            data(i,:) = channel_data; %#ok<AGROW>
        end

        % Call the single channel spectrogram function
        SleepEEG_singleChanSpec(subID, fnsuffix, data, commandDict.originFs, commandDict.project, srate, false, fileID, outputDir);
    end

    %%
    if commandDict.verbose
        disp([newline, mHead, 'Total time taken for SleepEEG_call()...'])
        fprintf(mSpace)
        toc(procstart)
    end

end

diary off
close all

end
