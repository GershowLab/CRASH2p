function varargout = annotateBehavioralStateGUI(varargin)
% ANNOTATEBEHAVIORALSTATEGUI MATLAB code for annotateBehavioralStateGUI.fig
%      ANNOTATEBEHAVIORALSTATEGUI, by itself, creates a new ANNOTATEBEHAVIORALSTATEGUI or raises the existing
%      singleton*.
%
%      H = ANNOTATEBEHAVIORALSTATEGUI returns the handle to a new ANNOTATEBEHAVIORALSTATEGUI or the handle to
%      the existing singleton*.
%
%      ANNOTATEBEHAVIORALSTATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANNOTATEBEHAVIORALSTATEGUI.M with the given input arguments.
%
%      ANNOTATEBEHAVIORALSTATEGUI('Property','Value',...) creates a new ANNOTATEBEHAVIORALSTATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before annotateBehavioralStateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to annotateBehavioralStateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help annotateBehavioralStateGUI

% Last Modified by GUIDE v2.5 15-May-2018 16:13:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @annotateBehavioralStateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @annotateBehavioralStateGUI_OutputFcn, ...
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


% --- Executes just before annotateBehavioralStateGUI is made visible.
function annotateBehavioralStateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to annotateBehavioralStateGUI (see VARARGIN)

% Choose default command line output for annotateBehavioralStateGUI
handles.output = hObject;
handles.tmr = varargin{1};
handles = initHandles(hObject,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes annotateBehavioralStateGUI wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = annotateBehavioralStateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.tmr;
delete(hObject);


% --- Executes on slider movement.
function videoPosition_Callback(hObject, eventdata, handles)
% hObject    handle to videoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
hObject.Value = round(hObject.Value);
% if (hObject.Value == handles.currFrame + 2)
%     hObject.Value = handles.currFrame + 1; %don't allow random frame skip -- don't know if this is needed or nto
% end
if (hObject.Value < handles.currFrame) %turn off autoset on backwards moves
    handles.persist_tags.Value = false;
end

handles = goToFrame(hObject, handles, hObject.Value, true);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function videoPosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to videoPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in set_state.
function set_state_Callback(hObject, eventdata, handles)
% hObject    handle to set_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% handles.sstic = tic();
% handles = setBehavioralState(hObject, handles);
% handles = goToFrame(hObject, handles, min(handles.currFrame + 1, handles.nframes), false);
% guidata(hObject, handles);

% --- Executes on button press in persist_tags.
function persist_tags_Callback(hObject, eventdata, handles)
% hObject    handle to persist_tags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of persist_tags

function handles = initHandles(hObject, handles)
handles.currFrame = handles.tmr.video.startFrame;
handles.persistTags = false;
handles.nframes = length(handles.tmr.video.et); %size(handles.tmr.videocube,3);
if (~(isa(handles.tmr.bsl, 'BehavioralStateList')) || isempty(handles.tmr.bsl.stateNames))
    handles.tmr = handles.tmr.initBSL;
end
if (~isfield(handles.tmr.video, 'micronDisplacement'))
    stagelpfreq = 5;
    handles.tmr.video.micronDisplacement = interp1(handles.tmr.tx, lowpass1D(handles.tmr.tracker.stageloc(1:2,:), 1000/(stagelpfreq*7.3))', handles.tmr.video.tx)';
end
if (isempty(handles.tmr.videocube))
    try
        handles.tmr = handles.tmr.openVideo();
    catch
        error ('no video loaded & cannot find video file');
    end
end
f = handles.tmr.getVidFrame(handles.tmr.video.startFrame);
handles.xl = [1 size(f, 2)];
handles.yl = [1 size(f, 1)];
handles.xl_full = handles.xl;
handles.yl_full = handles.yl;
handles.videoPosition.Min = 1;
handles.videoPosition.Max =handles.nframes;
handles.videoPosition.SliderStep = 1/(handles.nframes - 1)*[1,10];
handles = setBehavioralLabels(hObject, handles);
handles = goToFrame(hObject, handles, handles.currFrame, true);
guidata(hObject, handles);

function handles = setBehavioralLabels (hObject, handles)
if (isfield(handles, 'pushbuttons'))
    for j = 1:length(handles.pushbuttons)
        if (~ischar(handles.pushbuttons(j)))
            delete(handles.pushbuttons(j));
        end
    end
    handles.pushbuttons = [];
end
sn = handles.tmr.bsl.stateNames;
handles.behavioralLabels.Units = 'pixels';
ip = handles.behavioralLabels.Position;
w = ip(3);
h = ip(4)/(length(sn) + 1);
h0 = h/2;
for j = 1:length(sn)
    handles.pushbuttons(j) = uicontrol('Style', 'checkbox', 'String', sn{j}, 'Parent', handles.behavioralLabels, 'Position', [0 h0+(length(sn) - j)*h w h], 'Callback', @behavioralButtonCallback);
end
guidata(hObject,handles);

function handles = goToFrame(hObject, handles, newFrameNum, resetpushbuttons)
if (nargin < 3)
    newFrameNum = handles.currFrame;
end
if (nargin < 4)
    resetpushbuttons = false;
end

if (handles.persist_tags.Value)
    handles = setBehavioralState(hObject, handles,newFrameNum);
    resetpushbuttons = false;
end

handles.currFrame = newFrameNum;
handles.frame_number_text.String = num2str(handles.currFrame);
set(handles.videoPosition, 'Value', newFrameNum);
if (resetpushbuttons)
    handles = setBehavioralPushbuttons(hObject, handles, newFrameNum);
end
handles.tmr.drawFrame(newFrameNum, handles.videoAxes, 'overlayspot', false);
set(handles.videoAxes, 'XLim', handles.xl, 'YLim', handles.yl);
caxis(handles.videoAxes, round([handles.minZ.Value handles.maxZ.Value]));
axis(handles.videoAxes, 'equal');
handles.xl = handles.videoAxes.XLim;
handles.yl = handles.videoAxes.YLim;
handles = overlayText(hObject, handles);
guidata(hObject, handles);


function handles = overlayText(hObject, handles) 
s = get(handles.pushbuttons, 'String');
v = cell2mat(get(handles.pushbuttons, 'Value'));
s = s(logical(v));
try
    delete(handles.txt);
catch 
end
if (~isempty(s))
    s = regexprep(s, '_', ' ');
    handles.txt = text(handles.xl(2), handles.yl(2), s, 'Color', [0.2 1 1], 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 14, 'parent', handles.videoAxes);
end
guidata(hObject,handles);

function behavioralButtonCallback(hObject, eventData)
handles = guidata(hObject);
handles = overlayText(hObject,handles);
handles.persist_tags.Value = true;
guidata(hObject,handles);

function handles = setBehavioralState(hObject, handles, lastFrame)
s = get(handles.pushbuttons, 'String');
v = cell2mat(get(handles.pushbuttons, 'Value'));
s = s(logical(v));
existsAndDefault('lastFrame', handles.currFrame + 1);
for j = handles.currFrame:(lastFrame-1)
    handles.tmr.bsl = handles.tmr.bsl.setStateFrame(j, s);
end
autosaveTime = 600; %save every 10 mintues automatically
if (handles.autosave.Value)
    if (~isfield(handles, 'lastTime'))
        handles.lastTime = tic();
    end
    if (toc(handles.lastTime) > autosaveTime)
        handles.tmr.saveBSL();
        handles.lastTime = tic();
    end
end
guidata(hObject, handles);

function handles = setBehavioralPushbuttons(hObject, handles, newFrameNum)
for j = 1:length(handles.pushbuttons)
    set(handles.pushbuttons(j),'Value', handles.tmr.bsl.getStateFrame(newFrameNum, get(handles.pushbuttons(j), 'String')));
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function behavioralLabels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to behavioralLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
r = getrect(handles.videoAxes);
handles.xl = r(1) + [0 r(3)];
handles.yl = r(2) + [0 r(4)];
handles = goToFrame(hObject, handles);
guidata(hObject, handles);

% --- Executes on button press in unzoom.
function unzoom_Callback(hObject, eventdata, handles)
% hObject    handle to unzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.xl = handles.xl_full;
handles.yl = handles.yl_full;
handles = goToFrame(hObject, handles);
guidata(hObject, handles);


% --- Executes on slider movement.
function minZ_Callback(hObject, eventdata, handles)
% hObject    handle to minZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
goToFrame(hObject, handles);

% --- Executes during object creation, after setting all properties.
function minZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function maxZ_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to maxZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
goToFrame(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxZ_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to maxZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in autosave.
function autosave_Callback(hObject, eventdata, handles)
% hObject    handle to autosave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autosave


% --- Executes on button press in saveToDisk.
function saveToDisk_Callback(hObject, eventdata, handles)
% hObject    handle to saveToDisk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uiputfile(fullfile(fileparts(handles.tmr.video.filename), '*.csv'));
if (~isequal(filename, 0))
    handles.tmr.saveBSL(fullfile(pathname, filename));
    handles.lastTime = tic();
end
guidata(hObject, handles);

% --- Executes on button press in loadFromDisk.
function loadFromDisk_Callback(hObject, eventdata, handles)
% hObject    handle to loadFromDisk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tmr = handles.tmr.loadBSL();
guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, call UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over set_state.
function set_state_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to set_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isfield(handles, 'sstic'))
    handles.sstic = tic();
end
reptime = 0.05;
if (toc(handles.sstic) > reptime)
    set_state_Callback(hObject, eventdata, handles);
end

% --- Executes on key press with focus on set_state and none of its controls.
function set_state_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to set_state (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
