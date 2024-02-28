function [ dataDir, subject_list ] = SleepEEG_addpath(envpath)
%
% **ADSLEEPEEG_PREPROCESSING FUNCTION - ADDPATH**
%
% - used to facilitate addpath to various group code folders and Git repos
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Inputs:
%
%           - envpath:      a root path string returned by matlab system
%                           command matlabroot.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - dataDir:      the main directory with all subjects' data.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%%
if strcmp(envpath(1:13), '/Applications')
    % ---------- for Alex local env -------------
    % set up directories
    rootDir = '~/Dropbox/MIT_HST/NSRL_Purdon_Prerau/EEG';
    
    % Add path to Alex's code folder
    addpath([rootDir, '/code/sleepeeg_code'])
    
    % Add path to ANT_interface_code folder
    addpath([rootDir, '/code/ANT_interface_code'])
    
    % Add path to ANT importer
    addpath([rootDir, '/code/ANT_interface_code/ANTeepimport1.13'])
    
    % Add path to EEGLAB
    addpath([rootDir, '/code/ANT_interface_code/eeglab14_1_2b'])
    
    % Add path to EDF_deidentifier
    addpath([rootDir, '/code/sleepeeg_code/EDF_Deidentification_updated'])
    
    % Add path to Data folder
    currentpath = pwd;
    cd([rootDir '/data'])
    dataDir = pwd; %.mex doesn't accept ~/ in path
    addpath(dataDir)
    cd(currentpath)
    
elseif strcmp(envpath(1:7), '/autofs')
    % ---------- for NMR machine env -------------
    % for the ***adsleepeeg*** folder
    
    % set up root directories
    rootDir = '/autofs/cifs/adsleepeeg';
    
    % Add path to Alex's code folder
    addpath([rootDir, '/code/sleepeeg_code'])
    
    % Add path to ANT_interface_code folder
    addpath([rootDir, '/code/ANT_interface_code'])
    
    % Add path to ANT importer
    addpath([rootDir, '/code/ANT_interface_code/ANTeepimport1.13'])
    
    % Add path to EEGLAB
    addpath([rootDir, '/code/ANT_interface_code/eeglab14_1_2b'])
    
    % Add path to EDF_deidentifier
    addpath([rootDir, '/code/sleepeeg_code/EDF_Deidentification_updated'])
    
    % Add path to Data folder
    dataDir = [rootDir, '/archive/subject_data'];
    addpath(dataDir)
    
else
    % ---------- for ERIS cluster env -------------
    % for the ***adsleepeeg*** folder
    
    % set up root directories
    rootDir = '/data/adsleepeeg/';
    
    % Add path to Alex's code folder
    addpath([rootDir, '/code/sleepeeg_code'])
    
    % Add path to ANT_interface_code folder
    addpath([rootDir, '/code/ANT_interface_code'])
    
    % Add path to ANT importer
    addpath([rootDir, '/code/ANT_interface_code/ANTeepimport1.13'])
    
    % Add path to EEGLAB
    addpath([rootDir, '/code/ANT_interface_code/eeglab14_1_2b'])
    
    % Add path to EDF_deidentifier
    addpath([rootDir, '/code/sleepeeg_code/EDF_Deidentification_updated'])
    
    % Add path to Data folder
    dataDir = [rootDir, '/archive/subject_data'];
    addpath(dataDir)
    
end

if nargout > 1
    subject_list = dir(fullfile(dataDir, '*_*'));
    subject_list = {subject_list.name};
    subject_list = subject_list(cellfun(@(x) ~strcmp(x(1),'.'),subject_list));
    subject_list = subject_list';
end

end
