## SleepEEG_ Matlab functions
Author: Alex He 10/06/2020

  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Installation

In MATLAB, add path to the directory to this sleepeeg_code folder.

If on ERISOne cluster, SleepEEG_addpath.m will identify the adsleepeeg folder and add paths to the dependent functions such as ANT_interface_code folder.

Note: in order to use ANT import source codes on the Linux system, it requires GNU C-library version 2.14 or above. To use SleepEEG_call.m on ERISOne cluster, first ssh to one of nodes (yocto001, yocto002, grx03, grx04) and type in ldd --version to confirm glibc version. On these nodes, this function can be directly called from command window to perform various preprocessing actions on the raw .cnt data files.



## Usage

Once on a node with appropriate glibc. Type in terminal module load MATLAB/2019b (or above) to load Matlab to the session. Then Type matlab (use -nodesktop flag to ensure terminal session) to open a terminal MATLAB session. Add path to the location where the function SleepEEG_addpath.m is stored. In most cases, you want to use SleepEEG_call() as the wrapper function to call the other functions. Type  help SleepEEG_call for documentation on pre-defined settings.



## Function headers

=============================================================================

 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - ANTPOSSPEC**
 
  - used to plot the time traces of an anterior lead (84-Z2) and a
  posterior lead (25-Z10) and the multi-taper spectrograms. Outputs are
  saved to the same OutputDir
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - EEG:          data structure containing the EEG data and
                            impedance values.
 
            - chanlist:     Channel numbers to compute the spectrograms for.
                            default: [25, 84]
 
            - fileID:       name of the .cnt file (no .cnt suffix).
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - ADDPATH**
 
  - used to facilitate addpath to various group code folders and Git repos
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Inputs:
 
            - envpath:      a root path string returned by matlab system
                            command matlabroot.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - dataDir:      the main directory with all subjects' data.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - ANNOT**
 
  - this is a function used to organize the annotation structure extracted
  from Natus clinical sleep EEG system using the SleepEEG_loadedf()
  function
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Inputs:
 
            - header:       the header extracted by SleepEEG_loadedf().
 
            - annotation:   the annotation structure extracted by
                            SleepEEG_loadedf() but entries are in string
                            format so inconvenient to use.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - new_annotation:
                            an annotation structure containing the offset
                            duration, and annotation comments, with fields
                            cleaned up and organized.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - BRIDGE**
 
  - checks bridged electrodes during the EEG recording
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - EEG:          data structure containing the EEG data and
                            impedance values.
 
            - fileID:       name of the .cnt file (no .cnt suffix).
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - EB:           a structure containing various bridging
                            analysis outputs. This is also saved into the
                            analysis folder with the _bridging_EB_structure
                            suffix.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Command window display settings
  Beginning of command window messages.

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - BUILDANNOT**
 
  - this is a function used to construct an annotation channel signal based
  on an EDF+ file exported from the Natus clinical SleepWorks system at
  the MGH Sleep Center. The purpose of this is to add a dummy comment every
  5s so that they show up as annotations when the EDF is exported back
  into the Natus clinical SleepWorks system for us to extract the sleep
  stage scoring after the sleep tech has scored the imported study record
 
  This function automatically detects the number of samples in a single
  data record for the annotation channel and keeps it at that number. It
  will fill in the data records that don't have a comment in the original
  exported EDF+ file with dummy staging comments. As a result, The
  intervals between these dummy comments may not be exactly 5s, and will
  be shifted by the duration of a data record in order to maintain the same
  maximal samples in a single data record for the annotation channel.
 
  The output is a vector of signed integers of double type read as int16.
  You can treat this as the annotation channel data and feed it into
  blockEdfWrite() directly along with other channels that are also read as
  int16 (which is the default of SleepEEG_loadedf()).
 
  ===> Staging Comment abbreviated symbol inserted:

 
 
=============================================================================
 
 

 &#x1F536;
  **ADSLEEPEEG_PREPROCESSING FUNCTION - CALL**
 &#x1F536;
 
  - this is the master function to be called in the terminal window. One can
  call this function with subject ID and file suffix with various processing
  arguments that will be applied in an appropriate order
 
  Usage: [ EEG1 ] = <strong>SleepEEG_call</strong>(subID, fnsuffix, '<flag#1>',<arg#1>...'<flag#n>',<arg#n>);
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Declared Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep". if specified as
                            [] empty, then all .cnt files under the subID
                            directory will be processed.
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Optional Inputs:
 
            - []:           When <strong>SleepEEG_call</strong> is called without varargin
                            nor fnsuffix, only 'report' is true, all other
                            commands are set to false.
 
            - 'report':     To report ALL DATA files collected for the
                            subject, and prints in command window the types
                            of preprocessing that have been completed for
                            the subject's .cnt files.
                            default: false
 
            - 'verbose':    whether print messages during processing.
                            default: true
 
            - 'downsample': downsample to the a specified rate in Hz.
                            default: [true, 500 (Hz)]
 
            - 'dead':       whether to detect dead electrodes.
                            default: true
 
            - 'interp':     whether to interpolate dead electrodes.
                            default: false
 
            - 'reref':      whether re-reference to common average reference.
                            options are:
                                - AR (common average)
                                - Z3 (recording reference at Z3)
                                - [CH] (arbitrary channel number as reference)
                                - left mastoid
                                - linked mastoid
                                - REST (reference electrode standardization
                                        technique)
                                - contral mastoid
                                - LP (Laplacian reference based on duke layout)
                            default: {false, 'AR'}
 
            - 'oversave':   whether to overwrite existing .set files and
                            save the EEG structure again onto disk.
                            default: false
 
            - 'impedance':  whether check impedance values at the beginning
                            and the end of EEG recording.
                            default: true
 
            - 'bridge':     whether to do bridging analysis.
                            default: true
 
            - 'Pick2Spec':  channel number to plot the time traces,
                            locations, and spectrograms of 1 anterior and 1
                            posterior electrodes.
                            default: [25, 84]
 
            - 'originFs':   channel number to do spectrogram for two
                            electrodes in original sampling rate without
                            downsampling.
                            default: []
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - CONFIGDIR**
 
  - used to configure directories and specify path variables. Since lots of
  SleepEEG functions require to specify the right directories and places to
  find data, these actions are consolidated in this function
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep". if specified as
                            [] empty, then all .cnt files under the subID
                            directory will be returned as a cell array for
                            "fileID" output.
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - verbose:      whether print messages during processing.
                            default: true
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Outputs:
 
            - dataDir:      directory path containing all subjects' data.
 
            - datafn:       full path to the .cnt file including preceding
                            directories.
 
            - fileID:       name of the .cnt file (no .cnt suffix). if
                            fnsuffix input is specified as [] empty, then
                            all .cnt files under the subID directory will
                            be returned in datafn and subjectId as a cell
                            array.
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - CONSTANTCHAN**
 
  - Used to zero-center channels recording mostly constant values
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - EEG:          data structure containing the EEG data and
                            impedance values.
 
            - cthresh:      percentage threshold for identifying as a
                            constant channel.
                            default: 0.95
 
            - setzero:      whether to set the constant channels to zero.
                            default: false
 
            - ChanNum:      channel numbers to check for constant values.
                            default: [1]
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - Constantidx:  a vector of double type containing the channel
                            number of electrodes detected to have constant
                            values throughout the recording.
 
            - offsetlist:   a vector of double type with the same size as
                            the first input Constantidx containing the
                            offset values of the constant channels from
                            zero.
 
            - EEG:          an updated EEG structure if setzero was true
                            then EEG.data for the constant channels are set
                            to be be a constant zero recording.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - CONVERTSCORING**
 
  - this is a function used to read in the exported scoring text file after
  importing an EDF/EDF+ study back to the Natus system. After scoring is
  done by the sleep techs, the entire comment section is exported as a .txt
  file, and each comment is labelled with an epoch number and a scored
  sleep stage. We will read in the text file and update the Time column to
  be off-sets from the beginning of HD-EEG recordings. This converted
  scoring text file will then be imported with the HD-EEG data in BIDS
  format and stored as a MNE raw object. Various quality checks are
  completed to ensure accuracy of the extracted scoring
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - txtFN:        filename of the scoring .txt file.
 
            - edfFN:        filename of the exported .edf file that ends
                            with _deidentified.edf.
 
            - alignedfFN:   filename of the aligned .edf file that ends
                            with _aligned.edf.
 
            - filepath:     path to the "clinical" folder under the
                            subject's folder.
 
            - rvsFN:        filename of the reverse alignment .mat file
                            ouputted from trigger alignment code.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - DEAD**
 
  - Used to detect dead electrodes and electrodes with high impedances
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - EEG:          data structure containing the EEG data and
                            impedance values.
 
            - fileID:       name of the .cnt file (no .cnt suffix).
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
            - visualize:    whether to plot a topoplot of dead electrodes
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - Deadidx:      a vector of double type containing the channel
                            number of electrodes detected to be dead. The
                            reference electrode is excluded from this list
                            because it has constant zero recording by
                            definition.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - DOWNSAMPLE**
 
  - loads a recording from a .cnt file and downsamples the eego mylab high
  density EEG recording to 'dsrate' Hz in preparation for alignment with
  clinical PSG recordings
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - dsrate:       whether to downsample and desired downsampling rate.
                            default: [true, 500]
 
            - verbose:      whether print messages during processing.
                            default: true
 
            - oversave:     whether to overwrite existing .set files and
                            save the EEG structure again onto disk
                            default: false
 
            - dataDir:      directory path containing all subjects' data.
 
            - datafn:       full path to the .cnt file including preceding
                            directories. Can be obtained from
                            SleepEEG_configDir(subID, fnsuffix, verbose)
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - EEG:          an EEGLAB structure containing all information
                            of the recording in .cnt file.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - EXTRACTCHAN**
 
  - function used to extract several channels from overnight sleep EEG
  data recorded with ANT systems (asalab or eego lab)
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - channelNum:   a vector containing channel numbers to extract,
                            e.g. [25, 84].
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Outputs:
 
            - channelNum:   a vector containing channel numbers to extract,
                            repeating so that reconChan can take as input.
 
            - srate:        sampling rate of the extracted EEG channel, of
                            double type
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - FASTSCAN**
 
  - function used to prepare the digitization of electrode locations from
  Polhemus FastScanII scannner into data arrays and electrode labels to be
  read by mne.channels.read_dig_montage to create the montage for
  assembling a _raw.fif file for input in mne.gui.coregistration when
  generating the -trans.fif file during forward modeling
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            fastscan file, usually "night1(2)_fastscan".
 
                            *** the fastscan file name is
                            subID_fnsuffix.txt in the diretory '/fastscan'***
 
            - dataDir:      directory path containing all subjects' data.
                            default: will use SleepEEG_configDir()
 
            - criteria:     a 1x2 vector specifying the tolerable ranges of
                            angle deviation and distance deviation in
                            transforming from FastScan native coordinate
                            system to the Duke coordinate system.
                            defaul: [2, 0.2]
 
            - verbose:      whether print messages and plotting during
                            processing subject's electrode location files
                            default: true
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - fastscan:     a structure containing all information of the
                            digitization of electrodes acquired by Polhemus
                            FastScanII scanner. The fields of this
                            structure are as following:
 
                            fastscan.head
                                - point cloud of head surface
                            fastscan.electrode
                                - coordinates of electrodes in native
                                FastScan coordinate system
                            fastscan.landmark
                                - coordinates of right preauricular point,
                                nasion, and left preauricular point in
                                native FastScan coordinate system
                            fastscan.electrode_dukexyz
                                - coordinates of electrodes in the new Duke
                                Waveguard configuration coordinate system
                            fastscan.landmark_dukexyz
                                - coordinates of right preauricular
                                point, nasion, and left preauricular point
                                in the new Duke Waveguard configuration
                                coordinate system
                            fastscan.fschanlocs
                                - a structure containing the electrode
                                coordinates in the new Duke Waveguard
                                configuration coordinate system to be
                                integrated with an EEGLAB structure (EEG)
                                for the field chanlocs.
                            fastscan.elc_labels
                                - names of the channels in both
                                fastscan.electrode and
                                fastscan.electrode_dukexyz
                            fastscan.landmark_labels
                                - names of the fiducial landmarks in both
                                fastscan.landmark and
                                fastscan.landmark_dukexyz
                            fastscan.chanlocs_duke
                                - a structure containing the electrode
                                coordinates of the Duke Waveguard template
                                contained in the field chanlocs in an
                                EEGLAB structure (EEG) obtained from
                                calling ANT_interface_readcnt.m function.
                            fastscan.chanlocs_duke_reord
                                - a structure containing the electrode
                                coordinates of the Duke Waveguard template
                                contained in the field chanlocs in an
                                EEGLAB structure (EEG) obtained from
                                calling ANT_interface_readcnt.m function,
                                re-ordered to the same ordering as the
                                order of electrodes marked in Polhemus
                                FastScanII software (the order of
                                fastscan.electrode,
                                fastscan.electrode_dukexyz, and
                                fastscan.elc_labels).
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING SCRIPT - GENDOC**
 
  - this is the function used to generate documentations of all other
  functions. Not used for processing

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - HYPNO**
 
  - this is a function used to extract the scored sleep stages from
  exported .txt files from the clinical Natus system and plot hypnograms of
  various sorts
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            *** this function will use subID and fnsuffix
                                to load the .set file in the set folder if
                                'eeg' or 'both' is used as datasource:
                                subID_fnsuffix_ds500_Z3.set ***
 
            - channel:      a string character or a cell with different
                            channels when datasource is 'both'. If
                            datasource is 'eeg', then channel can be either
                            a string or a double type channel number. If
                            datasource is 'edf', then channel must be a
                            string character from F3,F4,C3,C4,O1,O2.
                            default: {'Z6', 'C3'}, or {87, 'C3'}
 
            - datasource:   datasource to plot the spectrogram from. Can be
                            'eeg' for grabbing channel from HD-EEG .set
                            file, 'edf' for grabbing channel from aligned
                            EDF+ file, or 'both' for plotting two
                            spectrograms, one for each datasource.
                            default: 'both'
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - IMPCHECK**
 
  - checks impedances at the beginning and the end of EEG recording
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - EEG:          data structure containing the EEG data and
                            impedance values.
 
            - fileID:       name of the .cnt file (no .cnt suffix).
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
            - visualize:    whether to generate plots of impedances
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - High_imp_electrode:
                            a vector of double type containing the channel
                            number of the electrodes with high impedance
                            values.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - LOADEDF**
  
  - this is a function used to load EDF/EDF+ file exported from the Natus
  clinical sleep EEG system at the MGH Sleep Center
  
  This function is based on the blockEdfLoad() function from Dr. Dennis A.
  Dean, II, Ph.D. Modifications are made in order to accommodate the 'EDF
  Annotations' channel exported from the clinical system. The rest of the
  header, signal header, and signal vector specfications of the exported
  EDF/EDF+ file are the same as the original EDF format (Kemp et al. 1992).
 
  However, the content largely follows the 12 specifications recommended by
  EDF+, although some of the channel fields like "transducer_type" and
  "prefiltering" are not reported.
  
  Check <a href="matlab:web https://www.edfplus.info/specs/edfplus.html">https://www.edfplus.info/specs/edfplus.html</a> for full specs.
 
  ##### Function Prototypes:
 
                                 header = <strong>SleepEEG_loadedf</strong>(edfFN)
                 [header, signalHeader] = <strong>SleepEEG_loadedf</strong>(edfFN)
     [header, signalHeader, signalCell] = <strong>SleepEEG_loadedf</strong>(edfFN)
     [header, signalHeader, signalCell] = <strong>SleepEEG_loadedf</strong>(edfFN, signalLabels)
     [header, signalHeader, signalCell] = <strong>SleepEEG_loadedf</strong>(edfFN, signalLabels, epochs)
 
  ##### Additional output (4th varargout):
 
          --- annotation   this is a structure containing the onset time,
                           duration, and content of various annotations
                           made by sleep techs during overnight recording.
 
  ##### To load just the annotations from the EDF+ file:
 
     [header, signalHeader, signalCell, annotation] = <strong>SleepEEG_loadedf</strong>(edfFN, {})
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - LOADSET**
 
  - Used to load .set EEG files based on subject ID and file suffix
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep". 
 
                            ***the final file name is in the form: subID_fnsuffix_ds500.set***
 
            - downsampled:  whether to add '_ds500' when loading data.
                            default: true
 
            - reref:        a string specifying the reference scheme of the
                            data to be loaded.
                            default: 'Z3'
 
            - verbose:      whether print messages during processing.
                            default: true
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - EEG:          an EEGLAB structure containing all information
                            of the recording in .set file.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - PLOT2DELC**
 
  - plots the 2D topoplot of EEG electrodes 
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - chanlocs:     a EEGLAB chanloc structure 
 
            - EEG:          data structure containing the EEG data and
                            chanloc.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - READSTAGING**
 
  - this is a function used to read in the exported staging data based on
  the manually introduced 'Staging Comment' annotations during importing an
  EDF/EDF+ study back to the Natus system. After scoring is done by the
  sleep techs, the entire comment section is exported as a .txt file, and
  each comment is labelled with an epoch number and a scored sleep stage.
  We will read in the text file and extract the sleep staging scoring as
  output. Various quality checks are completed to ensure accuracy of the
  extracted scoring
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - txtFN:        filename of the scoring .txt file.
 
            - edfFN:        filename of the exported .edf file that ends
                            with _deidentified.edf.
 
            - alignedfFN:   filename of the aligned .edf file that ends
                            with _aligned.edf.
 
            - filepath:     path to the "clinical" folder under the
                            subject's folder.
 
            - rvsFN:        filename of the reverse alignment .mat file
                            ouputted from trigger alignment code.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - stage_channel:
                            a vector containing the staging scoring in
                            the timeframe of the aligned EDF+ file.
                            I.e, it has the same number of samples as a
                            channel in the _aligned.edf file.
 
            - EEGvec_stage_channel:
                            a vector containing the staging scoring in the
                            timeframe of the EEG .set file. It has been
                            truncated from the stage_channel using the
                            reverse alignment structure outputted by the
                            trigger alignment code.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - RECONCHAN**
 
  - function used to extract one occipital channel from overnight sleep EEG
  data recorded with ANT systems
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - channelNum:   a vector containing channel numbers to extract,
                            e.g. [25, 84]. Required.
 
            - dataDir:      directory path containing all subjects' data.
 
            - datafn:       full path to the .cnt file including preceding
                            directories. Can be obtained from
                            SleepEEG_configDir(subID, fnsuffix, verbose)
 
            - fileID:       name of the .cnt file (no .cnt suffix). 
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - channelfn:    a cell array of filenames of saved channel
                            data.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - REPORT**
 
  - used to report the available .cnt files under the subID directory and
  the length of each EEG recording
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - SINGLECHANSPEC**
 
  - used to plot the time traces and spectrograms of several channels
  without downsampling. Typically such data are the outputs of functions
  SleepEEG_extractChan.m followed by SleepEEG_reconChan.m
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
 
            - fnsuffix:     a string of suffix to identify a specific
                            EEG .cnt file, usually "night1(2)_Sleep".
 
                            ***the final file name is in the form: subID_fnsuffix.cnt***
 
            - data:         a M x N matrix containing data of one or more
                            channels without downsampling (could be
                            downsampled as well).
 
            - channelNum:   a 1 x M vector containing channel numbers.
 
            - srate:        original sampling rate, double type.
                            default: 1000 (Hz)
 
            - hrscale:      whether to use hour scale on time axis, default
                            is false, which uses seconds.
                            default: false
 
            - fileID:       name of the .cnt file (no .cnt suffix).
 
            - outputDir:    directory path to folder for dumping all
                            analyses outputs.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================
 
 
 
  **ADSLEEPEEG_PREPROCESSING FUNCTION - VET**
 
  - used to vet all files that should be present in a subject folder and
  checks over all file naming conventions. Produce a subID_vet_info.txt
  file in the log folder that summarizes missing files
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ##### Inputs:
 
            - subID:        a string of subject identifier: e.g. "SP001".
                            this function will automatically search for
                            relevant files in the subject folder
                            Filenames to vet are hardcoded.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ###### Outputs:
 
            - no output for this function.
 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 
 
=============================================================================

 



## Contributing
Pull requests are welcomed. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.



## License
[MIT](https://choosealicense.com/licenses/mit/)
