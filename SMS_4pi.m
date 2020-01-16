function varargout = SMS_4pi(varargin)
% SMS_4PI MATLAB code for SMS_4pi.fig
%      SMS_4PI, by itself, creates a new SMS_4PI or raises the existing
%      singleton*.
%
%      H = SMS_4PI returns the handle to a new SMS_4PI or the handle to
%      the existing singleton*.
%
%      SMS_4PI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMS_4PI.M with the given input arguments.
%
%      SMS_4PI('Property','Value',...) creates a new SMS_4PI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SMS_4pi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SMS_4pi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SMS_4pi

% Last Modified by GUIDE v2.5 15-Apr-2017 20:58:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SMS_4pi_OpeningFcn, ...
    'gui_OutputFcn',  @SMS_4pi_OutputFcn, ...
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


% --- Executes just before SMS_4pi is made visible.
function SMS_4pi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SMS_4pi (see VARARGIN)
path = [pwd '\'];
set(handles.pathMainfolder, 'String', path); %set current working directory
set(handles.panelRecon,'Visible','off');
set(handles.selectGetpos,'value',1);
% set(handles.selectRecon,'value',1);
% Choose default command line output for SMS_4pi
addpath('sCMOS_calibration_files\');
scmos_cali_file='readout1_tmpcalibration.mat';

handles.scmos_cali_file = scmos_cali_file;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SMS_4pi wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SMS_4pi_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browseMainfolder.
function browseMainfolder_Callback(hObject, eventdata, handles)
% hObject    handle to browseMainfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cdir = get(handles.pathMainfolder, 'String'); %get current path

cdir(end) = [];
c = cdir(end); %get last character of path to file
while c ~= '\' %get path only
    cdir(end) = []; %delete filename character
    if length(cdir) == 0 %break if no path present
        cdir = [pwd '\'];
        break;
    end
    c = cdir(end);
end

if exist(cdir,'dir' )==7
    pathname = uigetdir(cdir,'Please select the main folder');
else
    cdir = pwd;
    c = cdir(end); %get last character of path to file
    while c ~= '\' %get path only
        cdir(end) = []; %delete filename character
        if isempty(cdir) %break if no path present
            break;
        end
        c = cdir(end); %get last character of path to file
    end
    pathname = uigetdir(cdir,'Please select the main folder');
end
if get(handles.selectGetpos,'value')==1
    filenum_threshold=str2num(get(handles.fileNoThresh,'string'));
    if pathname~=0
        set(handles.pathMainfolder,'string',[pathname '\']);
        folders=subdir(pathname);
        delin=[];
        kk=1;
        
        % count how many folders to perform analysis im
        for ff=1:1:size(folders,2)
            files=dir([folders{ff},'\*.dcimg']);
            if numel(files)<filenum_threshold;
                delin(kk)=ff;
                kk=kk+1;
            end
        end
        set(handles.programStatus,'string',['A total of ', num2str(size(folders,2)) ' sub folders were found (including intermediate root folders).'])
%         fprintf(['A total of ', num2str(size(folders,2)) ' sub folders were found (including intermediate root folders).\n']);
        folders(delin)=[];
        if isempty(folders) % check if it contains subfolder
            set(handles.pathSubfolder,'string','');
            return
        else
            set(handles.pathSubfolder,'string',folders);
            set(handles.pathSubfolder,'value',1:length(folders));
            
        end
        
    end
elseif get(handles.selectRecon,'value')==1
    if pathname~=0
        set(handles.pathMainfolder,'string',[pathname '\']);
        files=dir([pathname '\*tmpresult*.mat']);
        fileName={};
        for i =1:numel(files)
            fileName{i} = [pathname '\' files(i).name];
        end
        set(handles.pathSubfolder,'string',fileName);
        set(handles.pathSubfolder,'value',1:length(fileName));
    end
end



% --- Executes on selection change in pathSubfolder.
function pathSubfolder_Callback(hObject, eventdata, handles)
% hObject    handle to pathSubfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pathSubfolder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pathSubfolder
if get(handles.selectGetpos,'value')==1
    index_selected = get(handles.pathSubfolder,'value');
    folderlist = get(handles.pathSubfolder,'string');
    if class(folderlist) == 'char'
        temp{1}=folderlist;
        folderlist = temp;
    end
    folder_selected = folderlist(index_selected);
    if length(index_selected)==1&& exist(folder_selected{1},'dir' )==7
        
        files=dir([folder_selected{1} '\*.dcimg']);
        imageName = files(1).name;
        if ~isempty(strfind(imageName,'_642_'))
            channel = '642';
        elseif ~isempty(strfind(imageName,'_561_'))
            channel = '561';
        else
            channel = 'not 561 or 642'
        end
        set(handles.programStatus,'string',['Current Channel:' channel])
    elseif length(index_selected)>1
        set(handles.programStatus,'string','Multifolder Selected')
    end
elseif get(handles.selectRecon,'value')==1
    index_selected = get(handles.pathSubfolder,'value');
    filelist = get(handles.pathSubfolder,'string');
    file_selected = filelist(index_selected);
    if length(index_selected)==1
        currentfile=file_selected{1};
    I=find(currentfile=='\',1,'last');
    fileName=currentfile(I+1:end);
    if ~isempty(strfind(fileName,'_642_'))
            channel = '642';
        elseif ~isempty(strfind(fileName,'_561_'))
            channel = '561';
        else
            channel = 'not 561 or 642'
        end
        set(handles.programStatus,'string',['Current Channel:' channel])
    elseif length(index_selected)>1
        set(handles.programStatus,'string','Multifolder Selected')
    end
    
end

% --- Executes during object creation, after setting all properties.
function pathSubfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathSubfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton2.


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
posfile = get(handles.pathSubfolder, 'String'); %get path to image file
if class(posfile) == 'cell' %multiple files selected
    posfile = posfile{1}; %get first entry
end
c = posfile(end); %get last character of path to image file
while c ~= '\' %get path only
    posfile(end) = []; %delete filename character
    if length(posfile) == 0 %break if no path present
        break;
    end
    c = posfile(end); %get last charakter of path to image file
end
posfile = uigetdir(posfile, 'Save Position Directory'); %let user select output directory
if posfile ~= 0
    posfile = [posfile '\']; %add backslash
    set(handles.pathoutput, 'String', posfile); %set path
    clear posfile; %clear variable
else
    return
    endposfile = get(handles.pathSubfolder, 'String'); %get path to image file
    if class(posfile) == 'cell' %multiple files selected
        posfile = posfile{1}; %get first entry
    end
    c = posfile(end); %get last character of path to image file
    while c ~= '\' %get path only
        posfile(end) = []; %delete filename character
        if length(posfile) == 0 %break if no path present
            break;
        end
        c = posfile(end); %get last charakter of path to image file
    end
    posfile = uigetdir(posfile, 'Save Position Directory'); %let user select output directory
    if posfile ~= 0
        posfile = [posfile '\']; %add backslash
        set(handles.pathoutput, 'String', posfile); %set path
        clear posfile; %clear variable
    else
        return
    end
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function det_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to det_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of det_thresh as text
%        str2double(get(hObject,'String')) returns contents of det_thresh as a double


% --- Executes during object creation, after setting all properties.
function det_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to det_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function llrthresh1_Callback(hObject, eventdata, handles)
% hObject    handle to llrthresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of llrthresh1 as text
%        str2double(get(hObject,'String')) returns contents of llrthresh1 as a double


% --- Executes during object creation, after setting all properties.
function llrthresh1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to llrthresh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_Callback(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double


% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function loadTranform_Callback(hObject, eventdata, handles)
% hObject    handle to loadTranform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iPALMgetTranform(handles);


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to phaseshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reconstruct.
function reconstruct_Callback(hObject, eventdata, handles)
% hObject    handle to reconstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

iPALMreconstruction(handles);


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function center1_Callback(hObject, eventdata, handles)
% hObject    handle to center1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center1 as text
%        str2double(get(hObject,'String')) returns contents of center1 as a double


% --- Executes during object creation, after setting all properties.
function center1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function center2_Callback(hObject, eventdata, handles)
% hObject    handle to center2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center2 as text
%        str2double(get(hObject,'String')) returns contents of center2 as a double


% --- Executes during object creation, after setting all properties.
function center2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function center3_Callback(hObject, eventdata, handles)
% hObject    handle to center3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center3 as text
%        str2double(get(hObject,'String')) returns contents of center3 as a double


% --- Executes during object creation, after setting all properties.
function center3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function center4_Callback(hObject, eventdata, handles)
% hObject    handle to center4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center4 as text
%        str2double(get(hObject,'String')) returns contents of center4 as a double


% --- Executes during object creation, after setting all properties.
function center4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fileNoThresh_Callback(hObject, eventdata, handles)
% hObject    handle to fileNoThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileNoThresh as text
%        str2double(get(hObject,'String')) returns contents of fileNoThresh as a double


% --- Executes during object creation, after setting all properties.
function fileNoThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileNoThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function llrthresh2_Callback(hObject, eventdata, handles)
% hObject    handle to llrthresh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of llrthresh2 as text
%        str2double(get(hObject,'String')) returns contents of llrthresh2 as a double


% --- Executes during object creation, after setting all properties.
function llrthresh2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to llrthresh2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function photons_Callback(hObject, eventdata, handles)
% hObject    handle to photons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of photons as text
%        str2double(get(hObject,'String')) returns contents of photons as a double


% --- Executes during object creation, after setting all properties.
function photons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to photons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getPositions.
function getPositions_Callback(hObject, eventdata, handles)
% hObject    handle to getPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iPALMgetPosition(handles);


% --- Executes on button press in selectRecon.
function selectRecon_Callback(hObject, eventdata, handles)
% hObject    handle to selectRecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectRecon


% --- Executes on button press in selectGetpos.
function selectGetpos_Callback(hObject, eventdata, handles)
% hObject    handle to selectGetpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectGetpos


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(handles.uibuttongroup1.SelectedObject.String)



% --- Executes during object creation, after setting all properties.
function uipanel4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadPhaseshift_Callback(hObject, eventdata, handles)
% hObject    handle to loadPhaseshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iPALMgetPhaseShift(handles);


% --- Executes when uibuttongroup1 is resized.
function uibuttongroup1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function alignment_Callback(hObject, eventdata, handles)
% hObject    handle to alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in Selectionswitch.
function Selectionswitch_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Selectionswitch 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mainFolder = get(handles.pathMainfolder,'string');
if get(handles.selectGetpos,'value')==1
    set(handles.panelGetpos,'visible','on');
    set(handles.panelRecon,'visible','off');
    set(handles.programStatus,'string','');
    
    filenum_threshold=str2num(get(handles.fileNoThresh,'string'));
    if mainFolder(end)=='\'
        mainFolder(end)=[];
    end
    folders=subdir(mainFolder);
    delin=[];
    kk=1;
    
    % count how many folders to perform analysis im
    for ff=1:1:size(folders,2)
        files=dir([folders{ff},'\*.dcimg']);
        if numel(files)<filenum_threshold;
            delin(kk)=ff;
            kk=kk+1;
        end
    end
    fprintf(['A total of ', num2str(size(folders,2)) ' sub folders were found (including intermediate root folders).\n']);
    folders(delin)=[];
    if isempty(folders) % check if it contains subfolder
        set(handles.pathSubfolder,'string','');
        return
    else
        set(handles.pathSubfolder,'string',folders);
        set(handles.pathSubfolder,'value',1:length(folders));
        
    end
    
elseif get(handles.selectRecon,'value')==1
    position = get(handles.panelGetpos,'position');
    set(handles.panelRecon,'visible','on');
    set(handles.panelRecon,'position',position);
    set(handles.panelGetpos,'visible','off');
    set(handles.programStatus,'string','');
    
    
    files=rdir([mainFolder '**\*tmpresult*.mat']);
    fileName={};
    for i =1:numel(files)
        %         fileName{i} = [mainFolder files(i).name];
        fileName{i} = files(i).name;
    end
    set(handles.pathSubfolder,'string',fileName);
    set(handles.pathSubfolder,'value',1:length(fileName));
    
end



function CRLB_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to CRLB_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CRLB_thresh as text
%        str2double(get(hObject,'String')) returns contents of CRLB_thresh as a double


% --- Executes during object creation, after setting all properties.
function CRLB_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CRLB_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_thresh as text
%        str2double(get(hObject,'String')) returns contents of sigma_thresh as a double


% --- Executes during object creation, after setting all properties.
function sigma_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stopValue_Callback(hObject, eventdata, handles)
% hObject    handle to stopValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopValue as text
%        str2double(get(hObject,'String')) returns contents of stopValue as a double


% --- Executes during object creation, after setting all properties.
function stopValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mephiBin_Callback(hObject, eventdata, handles)
% hObject    handle to mephiBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mephiBin as text
%        str2double(get(hObject,'String')) returns contents of mephiBin as a double


% --- Executes during object creation, after setting all properties.
function mephiBin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mephiBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in phaseDriftCorr.
function phaseDriftCorr_Callback(hObject, eventdata, handles)
% hObject    handle to phaseDriftCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phaseDriftCorr



function numImages_Callback(hObject, eventdata, handles)
% hObject    handle to numImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numImages as text
%        str2double(get(hObject,'String')) returns contents of numImages as a double


% --- Executes during object creation, after setting all properties.
function numImages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in beadsDriftCorr.
function beadsDriftCorr_Callback(hObject, eventdata, handles)
% hObject    handle to beadsDriftCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beadsDriftCorr



function fit_flag_Callback(hObject, eventdata, handles)
% hObject    handle to fit_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_flag as text
%        str2double(get(hObject,'String')) returns contents of fit_flag as a double


% --- Executes during object creation, after setting all properties.
function fit_flag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MedFilter.
function MedFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MedFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MedFilter


% --- Executes on button press in Interpolation.
function Interpolation_Callback(hObject, eventdata, handles)
% hObject    handle to Interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Interpolation
