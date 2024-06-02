function varargout = EDF_Deidentifier(varargin)
% EDF_DEIDENTIFIER MATLAB code for EDF_Deidentifier.fig
%      EDF_DEIDENTIFIER, by itself, creates a new EDF_DEIDENTIFIER or raises the existing
%      singleton*.
%
%      H = EDF_DEIDENTIFIER returns the handle to a new EDF_DEIDENTIFIER or the handle to
%      the existing singleton*.
%
%      EDF_DEIDENTIFIER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDF_DEIDENTIFIER.M with the given input arguments.
%
%      EDF_DEIDENTIFIER('Property','Value',...) creates a new EDF_DEIDENTIFIER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EDF_Deidentifier_OpeningFcn gets called. An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EDF_Deidentifier_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EDF_Deidentifier

% Last Modified by GUIDE v2.5 27-Jun-2019 00:30:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EDF_Deidentifier_OpeningFcn, ...
    'gui_OutputFcn',  @EDF_Deidentifier_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EDF_Deidentifier is made visible.
function EDF_Deidentifier_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EDF_Deidentifier (see VARARGIN)

% Choose default command line output for EDF_Deidentifier
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EDF_Deidentifier wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EDF_Deidentifier_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select_file_pushbutton.
function select_file_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_file_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = hObject;
% set callback of this button to open a dialogue box to choose files from
handles=select_files(hObject, handles);

% change status
if ~isempty(handles.edf_filenames)
    set(handles.status_text,'String','File(s) loaded. Ready to de-identify selected files');
    set(handles.run_pushbutton,'Enable','on');
end

guidata(hObject, handles);


function [handles]= select_files(hObject, handles)
handles.output = hObject;
% function to let user select files to deidentify
[edf_filenames,edf_pathname] = uigetfile('*.edf','Select the EDF Files to deidentify','MultiSelect','on');

% change current directory to the selected files folder
cd(edf_pathname);

% save the filenames in the handle
handles.edf_filenames=edf_filenames;


guidata(hObject, handles);


% function deidentify_files_edf(edf_filenames,hObject, handles)
% handles.output = hObject;
% 
% % gather new EDF file suffix from the user
% prompt = {'Enter suffix to be appended to the deidentified file names:'};
% dlg_title = 'New EDF File Name Suffix';
% num_lines = 1;
% defaultans = {'_deidentifed'};
% suffix = inputdlg(prompt,dlg_title,num_lines,defaultans);
% 
% % go through all the files and convert them to edf
% if iscell(edf_filenames) % multiple files
%     for i=1:length(edf_filenames)
%         set(handles.status_text,'String',['De-identifying ', edf_filenames{i},'...']);drawnow();
%         edf_deidentify(edf_filenames{i},[edf_filenames{i}(1:end-4),suffix{1,1},'.edf']);
%     end
% else % single file
%     set(handles.status_text,'String',['De-identifying ', edf_filenames,'...']);drawnow();
%     edf_deidentify(edf_filenames,[edf_filenames(1:end-4),suffix{1,1},'.edf']);
% end
% set(handles.status_text,'String','De-identification done! Select new files to de-identify');
% set(handles.select_file_pushbutton,'Enable','on');
% guidata(hObject, handles)


% --- Executes on button press in rename_checkbox.
function rename_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to rename_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rename_checkbox
handles.output = hObject;

rename_check=get(hObject,'Value');

if rename_check
    set(handles.suffix_textbox,'Enable','on');
    set(handles.suffix_textbox,'BackGroundColor',[1 1 1]);
    set(handles.suffix_textbox,'String','_deidentified');
else
    set(handles.suffix_textbox,'Enable','inactive');
    set(handles.suffix_textbox,'BackGroundColor',[0.502 0.502 0.502]);
    set(handles.suffix_textbox,'String','');
end

guidata(hObject, handles);


function suffix_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to suffix_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of suffix_textbox as text
%        str2double(get(hObject,'String')) returns contents of suffix_textbox as a double
handles.output = hObject;

% Get the suffix. If empty, then save the deidentified with the same name.
suffix=get(hObject,'String');

% save it in the handles
handles.suffix=suffix;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function suffix_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to suffix_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.suffix_textbox = hObject;

guidata(hObject, handles);


% --- Executes on button press in edf_checkbox.
function edf_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to edf_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edf_checkbox
handles.output = hObject;

edf_save_check=get(hObject,'Value');
mat_save_check=get(handles.mat_checkbox,'Value');

if ~(edf_save_check || mat_save_check)
    set(hObject,'Value',1);
end

guidata(hObject, handles);


% --- Executes on button press in mat_checkbox.
function mat_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to mat_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mat_checkbox
handles.output = hObject;

mat_save_check=get(hObject,'Value');
edf_save_check=get(handles.edf_checkbox,'Value');

if ~(edf_save_check || mat_save_check)
    set(hObject,'Value',1);
end

guidata(hObject, handles);

% --- Executes on button press in run_pushbutton.
function run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = hObject;

set(handles.select_file_pushbutton,'Enable','off');
set(hObject,'Enable','off');

% Get the suffix. If empty, then save the deidentified with the same name.
suffix=get(handles.suffix_textbox,'String');

% get the edf and mat checkbox values
edf_save_check=get(handles.edf_checkbox,'Value');
mat_save_check=get(handles.mat_checkbox,'Value');

edf_filenames=handles.edf_filenames;

% check whether or not to overwrite file
% overwrite_option = get(get(handles.radiobutton_buttongroup,'SelectedObject'),'Tag');
overwrite = handles.overwrite.Value;


% save file as edf or mat or both
if edf_save_check || mat_save_check
    % go through all the files and convert them to edf
    if iscell(edf_filenames) % multiple files
        %         progressbar;
        for i=1:length(edf_filenames)
            set(handles.status_text,'String',['De-identifying ', edf_filenames{i},'...']);drawnow();
            
            % get name of new file
            new_filename=[edf_filenames{i}(1:end-4),suffix,'.edf'];
            
            
            % overwrite the edf file with the deidentified data and then
            % rename the file if rename is selected
            if (overwrite)&&(~isempty(suffix))
                edf_deidentify(edf_filenames{i},edf_filenames{i});
                movefile(edf_filenames{i},new_filename,'f');
            else 
                edf_deidentify(edf_filenames{i},new_filename);
            end
            % check if mat file also needs to be generated
            if mat_save_check
                [data,hdr]=readedf(new_filename);
                save([new_filename(1:(end-4)),'.mat'],'data','hdr');
            end
            %             progressbar(i/length(filenames));
        end
    else % single file
        set(handles.status_text,'String',['De-identifying ', edf_filenames,'...']);drawnow();
        
        % get name of new file
        new_filename=[edf_filenames(1:end-4),suffix,'.edf'];
        
        % overwrite the edf file with the deidentified data and then
        % rename the file if rename is selected
        if (overwrite)&&(~isempty(suffix))
            edf_deidentify(edf_filenames,edf_filenames);
            movefile(edf_filenames,new_filename,'f');
        else
            edf_deidentify(edf_filenames,new_filename);
        end
        
        % check if mat file also needs to be generated
        if mat_save_check
            [data,hdr]=readedf(new_filename);
            save([new_filename(1:(end-4)),'.mat'],'data','hdr');
        end
    end
    
    % update status
    set(handles.status_text,'String','De-identification done! Select new files to de-identify');
    set(handles.select_file_pushbutton,'Enable','on');
    set(hObject,'Enable','off');
    
end

guidata(hObject, handles);


% --- Executes when selected object is changed in radiobutton_buttongroup.
function radiobutton_buttongroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in radiobutton_buttongroup
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.output = hObject;

% overwrite/not
overwrite_option=get(eventdata.NewValue,'Tag');

if strcmp(overwrite_option,'overwrite_edf_radiobutton')
    overwrite=1;
    set(handles.rename_checkbox,'Value',0);
    
else
    overwrite=0;
    set(handles.rename_checkbox,'Value',1);
end

% save selection in handles
handles.overwrite=overwrite;

% call the callback to rename so that when making a copy the the add suffix is
% automatically checked and set to _deidentified
rename_checkbox_Callback(handles.rename_checkbox, eventdata, handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function status_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.status_text = hObject;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function run_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.run_pushbutton = hObject;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function select_file_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_file_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


handles.select_file_pushbutton = hObject;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edf_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edf_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.edf_checkbox = hObject;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mat_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mat_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mat_checkbox = hObject;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function overwrite_edf_radiobutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overwrite_edf_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.overwrite = hObject;

guidata(hObject, handles);
