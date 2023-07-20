function [ fastscan ] = SleepEEG_fastscan(subID, fnsuffix, dataDir, criteria, verbose)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - FASTSCAN**
%
% - function used to prepare the digitization of electrode locations from
% Polhemus FastScanII scannner into data arrays and electrode labels to be
% read by mne.channels.read_dig_montage to create the montage for
% assembling a _raw.fif file for input in mne.gui.coregistration when
% generating the -trans.fif file during forward modeling
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ##### Inputs:
%
%           - subID:        a string of subject identifier: e.g. "SP001".
%
%           - fnsuffix:     a string of suffix to identify a specific
%                           fastscan file, usually "night1(2)_fastscan".
%
%                           *** the fastscan file name is
%                           subID_fnsuffix.txt in the diretory '/fastscan'***
%
%           - dataDir:      directory path containing all subjects' data.
%                           default: will use SleepEEG_configDir()
%
%           - criteria:     a 1x2 vector specifying the tolerable ranges of
%                           angle deviation and distance deviation in
%                           transforming from FastScan native coordinate
%                           system to the Duke coordinate system.
%                           defaul: [2, 0.2]
%
%           - verbose:      whether print messages and plotting during
%                           processing subject's electrode location files
%                           default: true
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - fastscan:     a structure containing all information of the
%                           digitization of electrodes acquired by Polhemus
%                           FastScanII scanner. The fields of this
%                           structure are as following:
%
%                           fastscan.head
%                               - point cloud of head surface
%                           fastscan.electrode
%                               - coordinates of electrodes in native
%                               FastScan coordinate system
%                           fastscan.landmark
%                               - coordinates of right preauricular point,
%                               nasion, and left preauricular point in
%                               native FastScan coordinate system
%                           fastscan.electrode_dukexyz
%                               - coordinates of electrodes in the new Duke
%                               Waveguard configuration coordinate system
%                           fastscan.landmark_dukexyz
%                               - coordinates of right preauricular
%                               point, nasion, and left preauricular point
%                               in the new Duke Waveguard configuration
%                               coordinate system
%                           fastscan.fschanlocs
%                               - a structure containing the electrode
%                               coordinates in the new Duke Waveguard
%                               configuration coordinate system to be
%                               integrated with an EEGLAB structure (EEG)
%                               for the field chanlocs.
%                           fastscan.elc_labels
%                               - names of the channels in both
%                               fastscan.electrode and
%                               fastscan.electrode_dukexyz
%                           fastscan.landmark_labels
%                               - names of the fiducial landmarks in both
%                               fastscan.landmark and
%                               fastscan.landmark_dukexyz
%                           fastscan.chanlocs_duke
%                               - a structure containing the electrode
%                               coordinates of the Duke Waveguard template
%                               contained in the field chanlocs in an
%                               EEGLAB structure (EEG) obtained from
%                               calling ANT_interface_readcnt.m function.
%                           fastscan.chanlocs_duke_reord
%                               - a structure containing the electrode
%                               coordinates of the Duke Waveguard template
%                               contained in the field chanlocs in an
%                               EEGLAB structure (EEG) obtained from
%                               calling ANT_interface_readcnt.m function,
%                               re-ordered to the same ordering as the
%                               order of electrodes marked in Polhemus
%                               FastScanII software (the order of
%                               fastscan.electrode,
%                               fastscan.electrode_dukexyz, and
%                               fastscan.elc_labels).
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 4
    criteria = [2,0.2];
    verbose = true;
elseif nargin < 5
    verbose = true;
end

if ~exist('dataDir', 'var') || isempty(dataDir)
    [dataDir, ~, fileID, ~] = SleepEEG_configDir(subID, fnsuffix, verbose);
else
    fileID = [subID, '_', fnsuffix];
end

%% Command window display settings
% Beginning of command window messages.
mHead = 'SleepEEG: ';
% Spaces that can be used to replace mHead for better alignment of messages.
mSpace = repmat(sprintf(' '), 1, length(mHead));

%%
disp('-----------------------------')
disp([mHead 'SleepEEG_fastscan()']);
disp('-----------------------------')

%% Configure filenames and filepath
fsfn = [fileID, '_dig.mat'];
pcfn = [fileID, '.mat'];
markfn = [fileID, '.txt'];
subfolder = 'fastscan';
filepath = fullfile(dataDir, subID, subfolder);

%% Process and vet the FastScan coordinates
fastscan = ANT_MNE_fastscan(pcfn, markfn, fsfn, filepath, criteria, verbose);

end
