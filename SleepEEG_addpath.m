function [ dataDir, subject_list ] = SleepEEG_addpath(envpath, project)
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
%           - project:      an optional string to specify project name for
%                           path configuration.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ###### Outputs:
%
%           - dataDir:      the main directory with all subjects' data.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 2
    project = '';
end

%%
if strcmp(envpath(1:13), '/Applications')
    % ---------- for Alex local env -------------
    % set up directories
    if isempty(project)
        rootDir = '~/Dropbox/Active_projects/EEG';
    else
        rootDir = ['~/Dropbox/Active_projects/', project];
    end
    code_rootDir = '~/Dropbox/Active_projects/EEG/code';

    % Add path to Alex's code folder
    addpath([code_rootDir, '/sleepeeg_code'])

    % Add path to ANT_interface_code folder
    addpath([code_rootDir, '/ANT_interface_code'])

    % Add path to ANT importer
    addpath([code_rootDir, '/ANT_interface_code/ANTeepimport1.13'])

    % Add path to EEGLAB
    addpath([code_rootDir, '/ANT_interface_code/eeglab14_1_2b'])
    eeglab nogui;

    % Add path to EEGLAB firfilt plugin for pop_resample()
    addpath([code_rootDir, '/ANT_interface_code/eeglab14_1_2b/plugins/firfilt1.6.2'])

    % Add path to helper functions
    addpath([code_rootDir, '/sleepeeg_code/helper_functions'])
    addpath([code_rootDir, '/sleepeeg_code/helper_functions/EDF_Deidentification_updated'])

    % Add path to Data folder
    currentpath = pwd;
    cd([rootDir '/data'])
    dataDir = pwd; %.mex doesn't accept ~/ in path
    addpath(dataDir)
    cd(currentpath)

else
    % ---------- for workstation cerulean -------------
    % for the *** aeas *** project folder

    % set up root directories
    if isempty(project)
        rootDir = '/remote/projects/aeas';
    else
        rootDir = ['/remote/projects/', project];
    end
    code_rootDir = '/remote/users/mh1/code';

    % Add path to Alex's code folder
    addpath([code_rootDir, '/sleepeeg_code'])

    % Add path to ANT_interface_code folder
    addpath([code_rootDir, '/ant_interface_code'])

    % Add path to ANT importer
    addpath([code_rootDir, '/ant_interface_code/ANTeepimport1.13'])

    % Add path to EEGLAB
    addpath([code_rootDir, '/ant_interface_code/eeglab14_1_2b'])
    eeglab nogui;

    % Add path to EEGLAB firfilt plugin for pop_resample()
    addpath([code_rootDir, '/ant_interface_code/eeglab14_1_2b/plugins/firfilt1.6.2'])

    % Add path to helper functions
    addpath([code_rootDir, '/sleepeeg_code/helper_functions'])
    addpath([code_rootDir, '/sleepeeg_code/helper_functions/EDF_Deidentification_updated'])

    % Add path to Data folder
    dataDir = fullfile(rootDir, 'subject_data');
    addpath(dataDir)

end

if nargout > 1
    subject_list = dir(fullfile(dataDir, '*_*'));
    subject_list = {subject_list.name};
    subject_list = subject_list(cellfun(@(x) ~strcmp(x(1),'.'),subject_list));
    subject_list = subject_list';
end

end
