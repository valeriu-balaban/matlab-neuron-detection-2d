function varargout = NeuronSelectionGUI(varargin)
% NEURONSLECTIONGUI MATLAB code for NeuronSlectionGUI.fig
%      NEURONSLECTIONGUI, by itself, creates a new NEURONSLECTIONGUI or raises the existing
%      singleton*.
%
%      H = NEURONSLECTIONGUI returns the handle to a new NEURONSLECTIONGUI or the handle to
%      the existing singleton*.
%
%      NEURONSLECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURONSLECTIONGUI.M with the given input arguments.
%
%      NEURONSLECTIONGUI('Property','Value',...) creates a new NEURONSLECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuronSlectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuronSlectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuronSlectionGUI

% Last Modified by GUIDE v2.5 24-Sep-2018 20:25:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuronSlectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuronSlectionGUI_OutputFcn, ...
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


function I_out = setupImage(I_in, P)
% Resets the least significant bit of the image to 0 excepts at the
% locations specified by P (col, row, radius) which are set to 1

    I_out = bitset(I_in, 1, 0);

    for k = 1:size(P, 1)
        X = P(k, 1);
        Y = P(k, 2);
        R = P(k, 3);
        
        col_start = max(1, X - R - 1);
        col_end   = min(size(I_in, 2), X + R + 1);
        
        row_start = max(1, Y - R - 1);
        row_end   = min(size(I_in, 1), Y + R + 1);
        
        I_out(row_start:row_end, col_start:col_end) = bitset(I_out(row_start:row_end, col_start:col_end), 1, 1);
    end
    
function RGB = displayNeurons(I, P, S)
% Highlights the regions of I (N x M) specified by P (K x 3, with X, Y, Radius) 
% using circles. S is the index of the selected neuron, 0 if none should be selected.
% Returns an image RGB (N x M X 3)

I_n = setupImage(I, P);

RGB = zeros(size(I_n, 1), size(I_n, 2), 3, 'uint8');
RGB(:, :, 1) = uint8(122 * bitget(I_n, 1) + I_n);
RGB(:, :, 2) = uint8(153 * bitget(I_n, 1) + I_n);
RGB(:, :, 3) = uint8(I_n);

% Color the selected region
if S >= 1 && S <= size(P, 1)
    X = P(S, 1);
    Y = P(S, 2);
    R = P(S, 3);

    col_start = max(1, X - R - 1);
    col_end   = min(size(I, 2), X + R + 1);

    row_start = max(1, Y - R - 1);
    row_end   = min(size(I, 1), Y + R + 1);

    RGB(row_start:row_end, col_start:col_end, 1) = I_n(row_start:row_end, col_start:col_end, 1) + 154;
    RGB(row_start:row_end, col_start:col_end, 2) = I_n(row_start:row_end, col_start:col_end, 1); 
    RGB(row_start:row_end, col_start:col_end, 3) = I_n(row_start:row_end, col_start:col_end, 1) + 122;
end


function RGB = highlightImage(I, P, S)
% Highlights the regions of I (N x M) specified by P (K x 3, with X, Y, Radius) 
% using circles. S is the index of the selected neuron, 0 if none should be selected.
% Returns an image RGB (N x M X 3)

    I_mask = zeros(size(I, 1) + 1000, size(I, 2) + 1000, 'uint8');
    I_sel = I_mask;

    for k = 1:size(P, 1)
        X = P(k, 1) + 500;
        Y = P(k, 2) + 500;
        R = P(k, 3);
        if k == S
            I_sel(Y - R : Y + R, X - R : X + R) = 1;
        else
            I_mask(Y - R : Y + R, X - R : X + R) = 1;
        end
    end

    I_mask = I_mask(501:end - 500, 501:end - 500);
    I_sel  = I_sel(501:end - 500, 501:end - 500);
    
    RGB = zeros(size(I, 1), size(I, 2), 3, 'uint8');
    
    RGB(:, :, 1) = uint8(122 * I_mask + 153 * I_sel + I);
    RGB(:, :, 2) = uint8(153 * I_mask + I);
    RGB(:, :, 3) = uint8(122 * I_sel + I);

    
% --- Executes just before NeuronSlectionGUI is made visible.
function NeuronSlectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuronSlectionGUI (see VARARGIN)

% Choose default command line output for NeuronSlectionGUI
handles.output = hObject;

I = zeros(550, 650);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using NeuronSlectionGUI.
if strcmp(get(hObject,'Visible'),'off')
    imshow(I, 'Parent', handles.Axes)
end

% UIWAIT makes NeuronSlectionGUI wait for user response (see UIRESUME)
% uiwait(handles.MainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = NeuronSlectionGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.MainWindow)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.MainWindow,'Name') '?'],...
                     ['Close ' get(handles.MainWindow,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.MainWindow)

% --- Executes on button press in Open.
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [file,path] = uigetfile({'*.tif; *.tiff'}, 'Select a tif stack file');
    
    if isequal(file,0)
        handles.ErrorMessage.String = 'No file selected';
    else        
        handles.Open.UserData = fullfile(path, file);
        
        % Get number of stacks and allocate the memory
        info          = imfinfo(fullfile(path, file));
        handles.Image = zeros(info(1).Height, info(1).Width, length(info), 'uint8');
        handles.ColoredImage = cell(1, length(info));
        handles.Table.UserData = cell(1, length(info));
        handles.TableLabel.UserData = [];
        
        % Load previously saved neurons if they existed
        [filepath,name, ~] = fileparts(handles.Open.UserData);
        if exist([fullfile(filepath, name), '.mat'], 'file') == 2
            load([fullfile(filepath, name), '.mat'], 'NeuronLocations');
            handles.Table.UserData = NeuronLocations;
        end
        
        for k=1:length(info)
            handles.ErrorMessage.String = sprintf('Loading image %i of %i...', k, length(info));
            drawnow;
            
            handles.Image(:,:, k) = imread(fullfile(path, file), k);
            
            % Color images based on the loaded neurons 
            if ~isempty(handles.Table.UserData{k})
                P = cell2mat(handles.Table.UserData{k});
                %handles.ColoredImage{k} = highlightImage(handles.Image(:, :, k), P, 0);
                handles.ColoredImage{k} = displayNeurons(handles.Image(:, :, k), P, 0);
            else
                handles.Table.UserData{k} = cell(0, 3);
                handles.ColoredImage{k} = handles.Image(:,:, k);
            end
        end
        
         % Setup the slider
        handles.ImageIndexSlider.Enable = 'on';
        handles.ImageIndexSlider.Min = 1;
        handles.ImageIndexSlider.Max = length(info);
        handles.ImageIndexSlider.Value = 1;
        handles.ImageIndexSlider.SliderStep = [1, 1] ./ (length(info) - 1);
        
        % Enable image sliders
        handles.VerticalSlider.Enable = 'on';
        handles.HorizontalSlider.Enable = 'on';
        
        ImageSlider_Callback(handles.ImageIndexSlider, eventdata, handles);
        
        handles.ErrorMessage.String = ['Image file ', file, ' loaded'];
    end
    
    % Update handles structure
    guidata(hObject, handles)



% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    
% --- Executes during object creation, after setting all properties.
function ImageIndexSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageIndexSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on scroll wheel click while the figure is in focus.
function MainWindow_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to MainWindow (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

newValue = handles.ImageIndexSlider.Value - eventdata.VerticalScrollCount;

if newValue >= handles.ImageIndexSlider.Min && newValue <= handles.ImageIndexSlider.Max
    handles.ImageIndexSlider.Value = newValue;
    
    % Update the image
    ImageSlider_Callback(handles.ImageIndexSlider, eventdata, handles)
end


% --- Executes on key release with focus on MainWindow or any of its controls.
function MainWindow_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to MainWindow (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)

disp([eventdata.Key, eventdata.Character, eventdata.Modifier])

if strcmp(eventdata.Key, 'leftarrow') || strcmp(eventdata.Key, 'rightarrow')
    
    if strcmp(eventdata.Key, 'leftarrow')
        newValue = handles.ImageIndexSlider.Value - 1;
    else
        newValue = handles.ImageIndexSlider.Value + 1;
    end
    
    if newValue >= handles.ImageIndexSlider.Min && newValue <= handles.ImageIndexSlider.Max
        handles.ImageIndexSlider.Value = newValue;

        % Update the image
        ImageSlider_Callback(handles.ImageIndexSlider, eventdata, handles)
    end
    
elseif strcmp(eventdata.Key, 'd')
    k = round(handles.ImageIndexSlider.Value);

    if ~isempty(handles.TableLabel.UserData) && handles.TableLabel.UserData(1) == k
        idx = handles.TableLabel.UserData(2);
        
        handles.TableLabel.UserData = [];
        
        handles.Table.Data(idx, :) = [];
        handles.Table.UserData{k} = handles.Table.Data;
        
        % Save User Data
        NeuronLocations = handles.Table.UserData;
        [filepath,name, ~] = fileparts(handles.Open.UserData);
        save([fullfile(filepath, name), '.mat'], 'NeuronLocations');
        
        P = cell2mat(handles.Table.Data);
        handles.ColoredImage{k} = displayNeurons(handles.Image(:, :, k), P, 0);
    
        ImageSlider_Callback(hObject, eventdata, handles);
    end
elseif strcmp(eventdata.Key, 'escape')
    k = round(handles.ImageIndexSlider.Value);
    
    if ~isempty(handles.TableLabel.UserData) && handles.TableLabel.UserData(1) == k
        P = cell2mat(handles.Table.Data);
        handles.TableLabel.UserData = [];
        handles.ColoredImage{k} = displayNeurons(handles.Image(:, :, k), P, 0);
        
        ImageSlider_Callback(hObject, eventdata, handles);
    end  
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function MainWindow_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to MainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~strcmp(handles.ImageIndexSlider.Enable, 'on')
    return
end

ref = handles.Axes.Position;

p = get(hObject, 'currentpoint');

if p(1) > ref(1) && p(1) <= ref(1) + ref(3) && p(2) > ref(2) && p(2) <= ref(2) + ref(4)
    hOffset = ref(3) * 2;
    vOffset = ref(4) * 2;

    imHeight = size(handles.Image, 1);
    imWidth  = size(handles.Image, 2);

    startRow = 1 + round((imHeight - vOffset - 1) * (1 - handles.VerticalSlider.Value));
    startCol = 1 + round((imWidth  - hOffset - 1) * handles.HorizontalSlider.Value);
    
    c = startCol + (p(1) - ref(1)) * 2;
    r = startRow + (ref(2) + ref(4) - p(2)) * 2;
    
    k = round(handles.ImageIndexSlider.Value);
    
    % Check if selecting a neuron
    P = cell2mat(handles.Table.Data);

    % Compute the index of the closest point
    if isempty(P)
        idx = 0;
    else
        [~, idx] = min(sqrt((P(:, 1) - c).^2 + (P(:, 2) - r).^2));
    end

    if idx && ((P(idx, 1) - c) <= P(idx, 3) || (P(idx, 2) - r) <= P(idx, 3))
        handles.TableLabel.UserData = [k, idx];
    else
        if ~isempty(handles.TableLabel.UserData) && handles.TableLabel.UserData(1) == k
            idx = handles.TableLabel.UserData(2);

            handles.Table.Data(idx, 1:2) = {c, r};
            handles.Table.UserData{k} = handles.Table.Data;
            
            handles.TableLabel.UserData = [];
        else
            handles.TableLabel.UserData = [];
            handles.Table.Data = [handles.Table.Data; {c, r, 20}];
            handles.Table.UserData{k} = handles.Table.Data; 
        end
    end
    
    
    % Save User Data
    NeuronLocations = handles.Table.UserData;
    [filepath,name, ~] = fileparts(handles.Open.UserData);
    save([fullfile(filepath, name), '.mat'], 'NeuronLocations');
    
    P = cell2mat(handles.Table.Data);
    
    % Update the colored images
    if ~isempty(handles.TableLabel.UserData) && handles.TableLabel.UserData(1) == k
        idx_selected = handles.TableLabel.UserData(2);
    else
        idx_selected = 0;
    end
    handles.ColoredImage{k} = displayNeurons(handles.Image(:, :, k), P, idx_selected);
    
    ImageSlider_Callback(hObject, eventdata, handles);
    
    % Update handles structure
    guidata(hObject, handles)
    
end


% --- Executes during object creation, after setting all properties.
function HorizontalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HorizontalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, ~, handles)
% hObject    handle to VerticalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if ~strcmp(handles.ImageIndexSlider.Enable, 'on')
    return
end

hOffset = handles.Axes.Position(3) * 2;
vOffset = handles.Axes.Position(4) * 2;

imHeight = size(handles.Image, 1);
imWidth  = size(handles.Image, 2);

startRow = 1 + round((imHeight - vOffset - 1) * (1 - handles.VerticalSlider.Value));
startCol = 1 + round((imWidth  - hOffset - 1) * handles.HorizontalSlider.Value);

k = round(handles.ImageIndexSlider.Value);

if ~isempty(handles.TableLabel.UserData) && handles.TableLabel.UserData(1) ~= k
    P = cell2mat(handles.Table.Data);
    k_old = handles.TableLabel.UserData(1);
    handles.TableLabel.UserData = [];
    handles.ColoredImage{k_old} = displayNeurons(handles.Image(:, :, k_old), P, 0);
end

% Update neurons position for the current stack
handles.Table.Data = handles.Table.UserData{k};

% Tif stack index should be an integer
handles.ImageIndexSlider.Value = k;
   
imshow(handles.ColoredImage{k}(startRow:startRow + vOffset, startCol:startCol + hOffset, :), ...
       'Parent', handles.Axes)

% imshow(handles.Image(startRow:startRow + vOffset, startCol:startCol + hOffset, ...
%        round(handles.ImageIndexSlider.Value)), 'Parent', handles.Axes)

handles.ImageLabel.String = sprintf('Image %i of %i', handles.ImageIndexSlider.Value, ...
                                    handles.ImageIndexSlider.Max);
                                
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function VerticalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VerticalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected cell(s) is changed in Table.
function Table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to Table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)



function RadiusText_Callback(hObject, eventdata, handles)
% hObject    handle to RadiusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RadiusText as text
%        str2double(get(hObject,'String')) returns contents of RadiusText as a double


% --- Executes during object creation, after setting all properties.
function RadiusText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RadiusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EdgeWidthText_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeWidthText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgeWidthText as text
%        str2double(get(hObject,'String')) returns contents of EdgeWidthText as a double


% --- Executes during object creation, after setting all properties.
function EdgeWidthText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgeWidthText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ThetaThresholdText_Callback(hObject, eventdata, handles)
% hObject    handle to ThetaThresholdText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThetaThresholdText as text
%        str2double(get(hObject,'String')) returns contents of ThetaThresholdText as a double


% --- Executes during object creation, after setting all properties.
function ThetaThresholdText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThetaThresholdText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Detect.
function Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

k = round(handles.ImageIndexSlider.Value);
I = handles.Image(:, :, k);
RADIUS_LIST         = eval(handles.RadiusText.String);
EDGE_WIDTH          = eval(handles.EdgeWidthText.String);
THETA_THRESHOLD     = eval(handles.ThetaThresholdText.String);

[position, radius] = FindNeurons(I , RADIUS_LIST, EDGE_WIDTH, THETA_THRESHOLD);

handles.Table.Data          = num2cell([position, radius]);
handles.Table.UserData{k}   = handles.Table.Data;

% Save User Data
NeuronLocations = handles.Table.UserData;
[filepath,name, ~] = fileparts(handles.Open.UserData);
save([fullfile(filepath, name), '.mat'], 'NeuronLocations');

P = cell2mat(handles.Table.Data);

% Update the colored images
handles.ColoredImage{k} = displayNeurons(handles.Image(:, :, k), P, 0);

ImageSlider_Callback(hObject, eventdata, handles);
