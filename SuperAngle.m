function varargout = SuperAngle(varargin)
% SUPERANGLE MATLAB code for SuperAngle.fig
%      SUPERANGLE, by itself, creates a new SUPERANGLE or raises the existing
%      singleton*.
%
%      H = SUPERANGLE returns the handle to a new SUPERANGLE or the handle to
%      the existing singleton*.
%
%      SUPERANGLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUPERANGLE.M with the given input arguments.
%
%      SUPERANGLE('Property','Value',...) creates a new SUPERANGLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SuperAngle_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SuperAngle_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SuperAngle

% Last Modified by GUIDE v2.5 21-Feb-2016 17:17:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SuperAngle_OpeningFcn, ...
                   'gui_OutputFcn',  @SuperAngle_OutputFcn, ...
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


% --- Executes just before SuperAngle is made visible.
function SuperAngle_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SuperAngle (see VARARGIN)

% Choose default command line output for SuperAngle
handles.output = hObject;
clc;

disp ('SuperAngle V1.20GUI for Super-X EDS')
disp ('-------------------------------------------------------------------')
disp ('This software calculates:')
disp ('(a) Effective solid angle under specimen tilt and shift condition')
disp ('(b) Counts ratio')
disp ('(c) Absolute counts')
disp ('(d) Composition analysis based on counts/ratio method.')
disp ('     ')
disp ('It is based on numeric approach taking effects from multiple')
disp ('detector geometry, Be holder absorption. holder frame shadowing,')
disp ('and specimen absorption into consideration.')
disp ('Support for single and multiple detectors configuration.')
disp ('For strong absorption system, please set sample inclination angle.')
disp ('Spurious X-ray is not modeled, but please keep this effect in mind.')
disp ('     ')
disp ('Modeled and Coded by Weizong Xu, LeBeau Group')
disp ('Please find details in "A Numerical Model for Multiple Detector')
disp ('Energy Dispersive X-ray Spectroscopy in the Transmission Electron')
disp ('Microscope", Ultramicroscopy, 2016, in press');
disp ('     ')
disp ('Email James M. LeBeau (jmlebeau@ncsu.edu) or') 
disp ('Weizong Xu (wxu4@ncsu.edu) for information and supports.')
disp ('-------------------------------------------------------------------')
% May 2015

%=======Coordinate========
%   Titan G2 ChemiSTEM
%           Y+
%            |
%       D1   |     D2  
%     135deg |   45deg
%           (Z)-----X+
%            
%       D4         D3
%     225deg     315deg
%            
%=========================
handles.dAngle = 0.2; %Accuracy for angular intergration, 0.5 is fairly good for speed, 0.2 is better for accuracy
set(handles.edit_detAngle, 'String', num2str(handles.dAngle));
[ handles.detector_para, handles.tot_Det_num ] = Detector_input( 'detector.xlsx' );
[ handles.angle_search, handles.detector_para] = Detector_setup( handles.detector_para, handles.dAngle );
update_detector_info(hObject, eventdata, handles)

handles.SpuriousX=zeros(10); %disable SpuriousX function

handles.t_chk = 1; % 1 --> constant thickness t, other value --> constant spot during tilt, i.e. t will change
set(handles.checkbox2,'Value',1);
filename_Sinput='specimen_startup.xlsx';
[ handles.sample_para ] = Specimen_setup( filename_Sinput, handles.t_chk );
%[ handles.sample_para ] = Specimen_setup( 'specimen_SrTiO3.xlsx', handles.t_chk );
if (handles.sample_para(12)~=1)
    set(handles.checkbox1,'Value',1);
    set(handles.edit2, 'String', num2str(handles.sample_para(17)));
end
set(handles.edit3, 'String', num2str(handles.sample_para(5)));

handles.chk_Shadow = 2;
set(handles.checkbox3,'Value',0);

set(handles.edit4, 'String', num2str(handles.sample_para(18)));
set(handles.edit5, 'String', num2str(handles.sample_para(9)));
set(handles.edit6, 'String', num2str(handles.sample_para(10)));
set(handles.edit7, 'String', num2str(handles.sample_para(21)));
set(handles.edit8, 'String', num2str(handles.sample_para(3)));
set(handles.edit9, 'String', num2str(handles.sample_para(4)));
set(handles.edit10, 'String', num2str(handles.sample_para(1)));
set(handles.edit11, 'String', num2str(handles.sample_para(2)));
set(handles.edit18, 'String', num2str(handles.sample_para(22)));
set(handles.edit19, 'String', num2str(handles.sample_para(23)));
set(handles.popupmenu1, 'Value', handles.sample_para(13)-2);
set(handles.popupmenu2, 'Value', handles.sample_para(15)-2);
set(handles.popupmenu3, 'Value', handles.sample_para(14));
set(handles.popupmenu4, 'Value', handles.sample_para(16));

%[ handles.holder_para, handles.holder_frame_para ] = Holder_setup( handles.chk_Shadow, 'holder_dapo_No1.xlsx', handles.sample_para );
[ handles.holder_para, handles.holder_frame_para ] = Holder_setup( handles.chk_Shadow, 'holder_FEI_LB.xlsx', handles.sample_para );
if (handles.holder_para(10)==1)
    set(handles.popupmenu5, 'Value', 1);
    set(handles.checkbox4,'Value',1);
else
    set(handles.popupmenu5, 'Value',2);
end

set(handles.edit12, 'String', num2str(handles.holder_para(1)));
set(handles.edit13, 'String', num2str(handles.holder_para(2)));
set(handles.edit14, 'String', num2str(handles.holder_para(4)));

if (handles.holder_para(7)==1)
    set(handles.checkbox5,'Value',1);
else
    set(handles.checkbox5,'Value',0);
end

if (handles.holder_para(8)==1)
    set(handles.popupmenu6, 'Value', 1);
else
    set(handles.popupmenu6, 'Value',2);
    set(handles.edit15, 'String', num2str(handles.holder_para(9)));
end
%Set tilt series along x-direction as default
set(handles.checkbox7,'Value',1)
set(handles.checkbox8,'Value',1)
set(handles.checkbox9,'Value',0)
handles.chkXY = 1;
handles.line_cal=1;

handles.exp_file='';
%handles.chkXY = -1;
handles.A_exp_counts=0;
handles.B_exp_counts=0;
handles.chk_print = 2;

if (strcmp (filename_Sinput, 'specimen_Ni3Al_demo.xlsx')) %demo mode, has Mo slot grid, setup here
    handles.holder_para(7)= 1; %grid_chk(Shadow by sample grid? 1-Yes 2-No)
    handles.holder_para(8)= 1; %type_grid(standard 1x2 slot? (1=yes others no))
    handles.holder_para(9)= 0; %open_diameter_grid(If others, input grid open diameter (in mm))
    set(handles.checkbox5,'Value',1);
    handles.exp_file='Exp_data_Ni3Al_demo_X.xlsx';
    disp('Demo Ni3Al specimen info is loaded, have fun!')
    uiwait(msgbox('Demo Ni3Al specimen is loaded, have fun!'));
end

handles.search_Deg=30;
set(handles.edit20, 'String', num2str(handles.search_Deg));
handles.d_Deg = 1;
set(handles.edit21, 'String', num2str(handles.d_Deg));
handles.search_Deg_2D=30;
handles.d_Deg_2D = 1;
handles.search_Range_2D=0.6; %in mm
set(handles.edit22, 'String', num2str(handles.search_Range_2D));
handles.d_Range_2D = 0.1;
set(handles.edit23, 'String', num2str(handles.d_Range_2D));

handles.chk_cal_line = 1;
handles.chk_cal_comp = 1;
handles.chk_cal_tilt2D = 1;
handles.chk_cal_shift2D =1;
handles.chk_display = 2;

[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SuperAngle wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SuperAngle_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function update_detector_info(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text_det_takeoff1, 'String', '');
set(handles.text_det_azimuth1, 'String', '');
set(handles.text_det_distance1, 'String', '');
set(handles.text_det_tilt1, 'String', '');
set(handles.text_det_diameter1, 'String', '');
set(handles.text_det_solid1, 'String', '');
set(handles.text_det_takeoff2, 'String', '');
set(handles.text_det_azimuth2, 'String', '');
set(handles.text_det_distance2, 'String', '');
set(handles.text_det_tilt2, 'String', '');
set(handles.text_det_diameter2, 'String', '');
set(handles.text_det_solid2, 'String', '');
set(handles.text_det_takeoff3, 'String', '');
set(handles.text_det_azimuth3, 'String', '');
set(handles.text_det_distance3, 'String', '');
set(handles.text_det_tilt3, 'String', '');
set(handles.text_det_diameter3, 'String', '');
set(handles.text_det_solid3, 'String', '');
set(handles.text_det_takeoff4, 'String', '');
set(handles.text_det_azimuth4, 'String', '');
set(handles.text_det_distance4, 'String', '');
set(handles.text_det_tilt4, 'String', '');
set(handles.text_det_diameter4, 'String', '');
set(handles.text_det_solid4, 'String', '');
set(handles.text_note, 'String', '');

if (handles.tot_Det_num>=1)
set(handles.text_det_takeoff1, 'String', num2str(handles.detector_para(1,4)));
set(handles.text_det_azimuth1, 'String', num2str(handles.detector_para(1,5)));
set(handles.text_det_distance1, 'String', num2str(handles.detector_para(1,3)));
set(handles.text_det_tilt1, 'String', num2str(handles.detector_para(1,6)));
set(handles.text_det_diameter1, 'String', num2str(handles.detector_para(1,7)));
set(handles.text_det_solid1, 'String', num2str(handles.detector_para(1,8)));
end

if (handles.tot_Det_num>=2)
set(handles.text_det_takeoff2, 'String', num2str(handles.detector_para(2,4)));
set(handles.text_det_azimuth2, 'String', num2str(handles.detector_para(2,5)));
set(handles.text_det_distance2, 'String', num2str(handles.detector_para(2,3)));
set(handles.text_det_tilt2, 'String', num2str(handles.detector_para(2,6)));
set(handles.text_det_diameter2, 'String', num2str(handles.detector_para(2,7)));
set(handles.text_det_solid2, 'String', num2str(handles.detector_para(2,8)));
end

if (handles.tot_Det_num>=3)
set(handles.text_det_takeoff3, 'String', num2str(handles.detector_para(3,4)));
set(handles.text_det_azimuth3, 'String', num2str(handles.detector_para(3,5)));
set(handles.text_det_distance3, 'String', num2str(handles.detector_para(3,3)));
set(handles.text_det_tilt3, 'String', num2str(handles.detector_para(3,6)));
set(handles.text_det_diameter3, 'String', num2str(handles.detector_para(3,7)));
set(handles.text_det_solid3, 'String', num2str(handles.detector_para(3,8)));
end

if (handles.tot_Det_num>=4)
set(handles.text_det_takeoff4, 'String', num2str(handles.detector_para(4,4)));
set(handles.text_det_azimuth4, 'String', num2str(handles.detector_para(4,5)));
set(handles.text_det_distance4, 'String', num2str(handles.detector_para(4,3)));
set(handles.text_det_tilt4, 'String', num2str(handles.detector_para(4,6)));
set(handles.text_det_diameter4, 'String', num2str(handles.detector_para(4,7)));
set(handles.text_det_solid4, 'String', num2str(handles.detector_para(4,8)));
end

if (handles.tot_Det_num>=5)
set(handles.text_note, 'String', 'More than 4 detectors');
end


function update_FEI_holder(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EleA_num = handles.sample_para(13);
EleA_shell = handles.sample_para(14); %1) K; 2) L; 3) M;
EleB_num = handles.sample_para(15);
EleB_shell = handles.sample_para(16); %1) K; 2) L; 3) M;
Absorp_table (:,:,1) = xlsread ('Absorption coefficient_K.xlsx'); %K-shell cm2/g
Absorp_table (:,:,2) = xlsread ('Absorption coefficient_L.xlsx'); %L-shell cm2/g
Absorp_table (:,:,3) = xlsread ('Absorption coefficient_M.xlsx'); %M-shell cm2/g
Be_density = 1.848;%g/cm3
uA_holder = Absorp_table(4,EleA_num,EleA_shell)*Be_density; %atomic number of Be is 4
uB_holder = Absorp_table(4,EleB_num,EleB_shell)*Be_density;
disp('Setup ideal holder parameters for FEI low background holder in Titan');

handles.holder_para(5)=uA_holder;
handles.holder_para(6)=uB_holder;
handles.holder_para(1)=0.22;
handles.holder_para(2)=3.10;
handles.holder_para(4)=20;
handles.chk_Shadow = 2;


set(handles.checkbox4,'Value',1);
set(handles.edit12, 'String', num2str(handles.holder_para(1)));
set(handles.edit13, 'String', num2str(handles.holder_para(2)));
set(handles.edit14, 'String', num2str(handles.holder_para(4)));
set(handles.checkbox3,'Value',0);
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)







function edit_detAngle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_detAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_detAngle as text
%        str2double(get(hObject,'String')) returns contents of edit_detAngle as a double
handles.dAngle=str2double(get(hObject,'String'));
disp('dAngle is set as ')
disp(handles.dAngle)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_detAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_detAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in input_detector.
function input_detector_Callback(hObject, eventdata, handles)
% hObject    handle to input_detector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.xlsx'},'Detector input file')
%[ handles.detector_para, handles.tot_Det_num ] = Detector_input( 'detector.xlsx' );
if (filename ~=0)
[ handles.detector_para, handles.tot_Det_num ] = Detector_input( filename );
[ handles.angle_search, handles.detector_para] = Detector_setup( handles.detector_para, handles.dAngle );
%set(handles.staticText1, 'String', num2str(value));
set(handles.file_det, 'String', filename);
update_detector_info(hObject, eventdata, handles);
guidata(hObject,handles)
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
[ elementA_number ] = periodic_info( str, val);
%[ elementA_weight ] = get_element_weight( elementA_number )
handles.sample_para(13)= elementA_number;
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
[ elementB_number ] = periodic_info( str, val);
%[ elementB_weight ] = get_element_weight( elementB_number )
handles.sample_para(15)= elementB_number;
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
% Save the handles structure.
guidata(hObject,handles)

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


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
case 'K' 
xray_shellA = 1;
case 'L' 
xray_shellA = 2;
case 'M' 
xray_shellA = 3;
case 'N' 
xray_shellA = 4;
end
handles.sample_para(14)= xray_shellA;
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
case 'K' 
xray_shellB = 1;
case 'L' 
xray_shellB = 2;
case 'M' 
xray_shellB = 3;
case 'N' 
xray_shellB = 4;
end
handles.sample_para(16)= xray_shellB;
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.sample_para(17)=str2double(get(hObject,'String'));
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
disp('Ideal atomic ratio of the specimen is ')
disp(handles.sample_para(17))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
chk = get(hObject,'Value');
if chk == 1
    handles.sample_para(12)=0; %cal_chk=1 direct cal from exp data, !=1) compare with ideal composition
    disp('calculation or compare using known composition')
else
    handles.sample_para(12)=1;
    disp('composition is unknown')
end
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
chk = get(hObject,'Value');
if chk == 1
    handles.t_chk = 1;
    disp('Same specimen thickness region')
else
    handles.t_chk = 0;
    disp('Same specimen location')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
chk = get(hObject,'Value');
if chk == 1
    handles.chk_Shadow = -1;
    handles.holder_para(3) = -1;
    disp('Set as no beam blocking. Danger! Not a real situation!')
else
    handles.chk_Shadow = 2;
    handles.holder_para(3) = 2;
    disp('Considering absorption by holders, real situation.')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.sample_para(5)=str2double(get(hObject,'String'));
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
disp('Specimen thickness is set as')
disp(handles.sample_para(5))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
handles.sample_para(18)=str2double(get(hObject,'String'));
[ handles.sample_para ] = update_specimen_para( handles.sample_para );
[ handles.holder_para ] = update_holder_para( handles.holder_para, handles.sample_para );
disp('Specimen density is set as')
disp(handles.sample_para(18))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
handles.sample_para(9)=str2double(get(hObject,'String'));
disp('Specimen shift in x-direction is set as')
disp(handles.sample_para(9))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
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
handles.sample_para(10)=str2double(get(hObject,'String'));
disp('Specimen shift in y-direction is set as')
disp(handles.sample_para(10))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
handles.sample_para(21)=str2double(get(hObject,'String'));
disp('Specimen depth in z-direction is set as')
disp(handles.sample_para(21))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.sample_para(3)=str2double(get(hObject,'String'));
disp('Specimen tilt in x-direction is set as')
disp(handles.sample_para(3))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
handles.sample_para(4)=str2double(get(hObject,'String'));
disp('Specimen tilt in y-direction is set as')
disp(handles.sample_para(4))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
handles.sample_para(1)=str2double(get(hObject,'String'));
disp('Incline angle of the specimen top surface with respect to holder x-axis is set as')
disp(handles.sample_para(1))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
handles.sample_para(2)=str2double(get(hObject,'String'));
disp('Incline angle of the specimen top surface with respect to holder y-axis is set as')
disp(handles.sample_para(2))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
EleA_num = handles.sample_para(13);
EleA_shell = handles.sample_para(14); %1) K; 2) L; 3) M;
EleB_num = handles.sample_para(15);
EleB_shell = handles.sample_para(16); %1) K; 2) L; 3) M;
Absorp_table (:,:,1) = xlsread ('Absorption coefficient_K.xlsx'); %K-shell cm2/g
Absorp_table (:,:,2) = xlsread ('Absorption coefficient_L.xlsx'); %L-shell cm2/g
Absorp_table (:,:,3) = xlsread ('Absorption coefficient_M.xlsx'); %M-shell cm2/g
chk = get(hObject,'Value');
if chk == 1
    Be_density = 1.848;%g/cm3
    uA_holder = Absorp_table(4,EleA_num,EleA_shell)*Be_density; %atomic number of Be is 4
    uB_holder = Absorp_table(4,EleB_num,EleB_shell)*Be_density;
    disp('Setup Be specimen carrier')
else
    Mo_density = 10.2;
    uA_holder = Absorp_table(42,EleA_num,EleA_shell)*Mo_density; %assume it is Mo (atom#42)
    uB_holder = Absorp_table(42,EleB_num,EleB_shell)*Mo_density;
    disp('No filtering effect from specimen carrier, X-ray will be fully blocked.')
end
handles.holder_para(5)=uA_holder;
handles.holder_para(6)=uB_holder;
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
case 'FEI low background holder' 
update_FEI_holder(hObject, eventdata, handles);
case 'Other' 
%xray_shellA = 2;
end
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
handles.holder_para(1)=str2double(get(hObject,'String'));
disp('Ideal sample depth in the specimen carrier is set as');disp(handles.holder_para(1));
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
handles.holder_para(2)=str2double(get(hObject,'String'));
disp('Open diameter of specimen carrier is set as')
disp(handles.holder_para(2))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
handles.holder_para(4)=str2double(get(hObject,'String'));
disp('The cone angle of the specimen carrier is set as')
disp(handles.holder_para(4))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
chk = get(hObject,'Value');
if chk == 1
    handles.holder_para(7) = 1;
    disp('X-ray can be shadowed by the supporting grid.')
else
    handles.holder_para(7) = 0;
    disp('X-ray will not be shadowed by the supporting grid.')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
case '1mmx2mm slot' 
handles.holder_para(8)=1;
disp('Setup grid as 1mmx2mm slot')
disp('Note: the length direction of the slot is set parallel to x-axis.') 
disp('Code modification is needed for other slot orientation setup.')
case 'Circular holes' 
handles.holder_para(8)=0;
disp('Setup grid as circular holes')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
handles.holder_para(9)=str2double(get(hObject,'String'));
disp('The inner diameter of a circular grid is set as');disp(handles.holder_para(9));
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

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


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
chk = get(hObject,'Value');
if chk == 1
    handles.chk_spot_cal = 1;
    disp('Enable calculation from single spot data.');
else
    handles.chk_spot_cal = 0;
    disp('Disable calculation from single spot data.');
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
chk = get(hObject,'Value');
if chk == 1
    set(handles.checkbox8,'Value',1)
    set(handles.checkbox9,'Value',0)
    handles.chkXY = 1;
    handles.line_cal=1;
else
    set(handles.checkbox8,'Value',0)
    set(handles.checkbox9,'Value',0)
    handles.line_cal=0;
    %set(handles.text_exp, 'String', 'Experiment data input');
    %handles.exp_file ='';
end
guidata(hObject,handles)

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8
chk = get(hObject,'Value');
if chk == 1 %check tiltX
    set(handles.checkbox9,'Value',0)
    set(handles.checkbox7,'Value',1)
    handles.chkXY=1;
    handles.line_cal=1;
else
    set(handles.checkbox9,'Value',1)
    handles.chkXY=2;
end
guidata(hObject,handles)

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
chk = get(hObject,'Value');
if chk == 1 %check tiltX
    set(handles.checkbox8,'Value',0)
    set(handles.checkbox7,'Value',1)
    handles.chkXY=2;
    handles.line_cal=1;
else
    set(handles.checkbox8,'Value',1)
    handles.chkXY=1;
end
guidata(hObject,handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*.xlsx'},'Experiment data input file')
if (filename ~=0)
set(handles.text67, 'String', filename);
handles.exp_file =filename;
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
handles.A_exp_counts=str2double(get(hObject,'String'));
disp('Experimental counts of element A is set as')
disp(handles.A_exp_counts)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
handles.B_exp_counts=str2double(get(hObject,'String'));
disp('Experimental counts of element B is set as')
disp(handles.B_exp_counts)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
handles.sample_para(22)=str2double(get(hObject,'String'));
disp('Probe current is set as')
disp(handles.sample_para(22))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
handles.sample_para(23)=str2double(get(hObject,'String'));
disp('Live acquisition time is set as')
disp(handles.sample_para(23))
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% call function single_spot( TiltX, TiltY, tot_Det_num, sample_para, holder_para, angle_search, SpuriousX)
% disp('Counts for single spot is under construction')
% [ cwText ] = GetCommandWindow();
% set(handles.edit24, 'String', cwText);
if (handles.sample_para(12)~=0)
    uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
    disp('Unable to calculate if composition is not known, please check specimen input file.');
else   
    tiltX=handles.sample_para(3);
    tiltY=handles.sample_para(4);
    disp('running...')
    uiwait(msgbox('Click OK to run'));
    msgbox('Running...');
    [ output_Point, output_omega ] = single_spot_counts(tiltX, tiltY, handles.tot_Det_num, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX); % counts*1e3;
    disp('Done.')

[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
handles.output_Point=output_Point;
handles.output_omega=output_omega;
end
guidata(hObject,handles)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.sample_para(12)~=0)
    uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
    disp('Unable to calculate if composition is not known, please check specimen input file.');
else    
    tiltX=handles.sample_para(3);
    tiltY=handles.sample_para(4);
    disp('running...')
    uiwait(msgbox('Click OK to run'));
    msgbox('Running...');
    [ output_Point, output_omega ] = single_spot(tiltX, tiltY, handles.tot_Det_num, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX); % counts*1e3;
    disp('Done.')

[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
handles.output_Point=output_Point;
handles.output_omega=output_omega;
end
guidata(hObject,handles)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.A_exp_counts >0 && handles.B_exp_counts >0)
    %disp(handles.A_exp_counts);disp(handles.B_exp_counts);
    disp('running...')
    uiwait(msgbox('Click OK to run'));
    msgbox('Running...');
    [handles.comp_ratio_weight_Spot, handles.comp_ratio_atomic_spot] = composition_cal_single_spot( handles.A_exp_counts, handles.B_exp_counts, handles.tot_Det_num, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
else
    uiwait(msgbox('Error! No experimental data input'));
    disp('Error! No experimental data input');
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.chkXY ==-1|| handles.sample_para(12)~=0)
    if (handles.chkXY ==-1)
        uiwait(msgbox('Unable to calculate if Tilt series direction is not set up.'));
        disp('Unable to calculate if Tilt series direction is not set up.');
    end
    if (handles.sample_para(12)~=0)
        uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
        disp('Unable to calculate if composition is not known, please check specimen input file.');
    end
else
    disp('running... It may take long time to finish...')
    uiwait(msgbox('Note: it may take long time to finish...'));
    msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[handles.Al_line, handles.Ni_line, handles.Al_Abso, handles.Ni_Abso, handles.Absrp_line ] = line_search( handles.chkXY, handles.tot_Det_num, handles.search_Deg, handles.d_Deg, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_line=2;
end
if (handles.chk_cal_line ==2)
line_display( handles.chk_print, handles.chkXY, handles.tot_Det_num, handles.exp_file, handles.search_Deg, handles.Al_line, handles.Ni_line, handles.Al_Abso, handles.Ni_Abso, handles.Absrp_line, handles.sample_para );
else
disp('Warning! No calculated data to display')
end
end
guidata(hObject,handles)
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.chkXY ==-1 || handles.sample_para(12)~=0)
    if (handles.chkXY ==-1)
        uiwait(msgbox('Unable to calculate if Tilt series direction is not set up.'));
        disp('Unable to calculate if Tilt series direction is not set up.');
    end
    if (handles.sample_para(12)~=0)
        uiwait(msgbox('Unable to calculate if composition is not known, please check specimen input file.'));
        disp('Unable to calculate if composition is not known, please check specimen input file.');
    end
else
    disp('running... It may take long time to finish...')
    uiwait(msgbox('Note: it may take long time to finish...'));
    msgbox('Running...');
if (handles.chk_display == 2)
[handles.Al_line, handles.Ni_line, handles.Al_Abso, handles.Ni_Abso, handles.Absrp_line ] = line_search( handles.chkXY, handles.tot_Det_num, handles.search_Deg, handles.d_Deg, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_line=2;
end
if (handles.chk_cal_line ==2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
line_display_Counts( handles.chk_print, handles.chkXY, handles.tot_Det_num, handles.exp_file, handles.search_Deg, handles.Al_line, handles.Ni_line, handles.Al_Abso, handles.Ni_Abso, handles.Absrp_line, handles.sample_para );
else
disp('Warning! No calculated data to display')
end
end
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.chkXY == -1 || isempty(handles.exp_file))
    if (handles.chkXY == -1)
        uiwait(msgbox('Error! Tilt series direction is not set up.'));
        disp('Error! Tilt series direction is not set up.')
    end
    
    if (isempty(handles.exp_file))
        uiwait(msgbox('Error! No experimental file input'));
        disp('Error! No experimental file input')
    end
else
    if (handles.sample_para(12)==0)
        uiwait(msgbox('Known composition! Deviation of prediction will be compared instead.'));
        disp('Known composition! Deviation of prediction will be compared instead.');
    end
disp('running... It may take long time to finish...')
uiwait(msgbox('Note: it may take long time to finish...'));
msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[handles.comp_ratio_out] = composition_cal( handles.exp_file, handles.tot_Det_num, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX, handles.chkXY );
    if (handles.sample_para(12)==0)
       [handles.comp_A_out, handles.comp_B_out] = composition_cal_absolute( handles.exp_file, handles.tot_Det_num, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX, handles.chkXY );
       handles.chk_cal_comp=2;
    else
       handles.chk_cal_comp=1; 
    end
end
if (handles.chk_cal_comp ==2)
[handles.diff_wt]=composition_display( handles.chk_print, handles.comp_ratio_out, handles.sample_para, handles.chkXY );
[handles.diff_wt_A, handles.diff_wt_B]=composition_display_absolute( handles.chk_print, handles.comp_A_out, handles.comp_B_out, handles.sample_para, handles.chkXY );
guidata(hObject,handles)
else
disp('No calculated data to display')
%uiwait(msgbox('No calculated data to display'));
end
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('running... It may take long time to finish...')
uiwait(msgbox('Note: it may take long time to finish...'));
msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[ handles.Al_map, handles.Ni_map, handles.ratio_map ] = XY_search( handles.tot_Det_num, handles.search_Deg_2D, handles.d_Deg_2D, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_tilt2D=2;
end
if (handles.chk_cal_tilt2D ==2)
XY_display2_solidAngle( handles.chk_print, handles.Al_map, handles.Ni_map, handles.ratio_map, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Deg_2D, handles.d_Deg_2D );
%XY_display2_counts( handles.chk_print, handles.Al_map, handles.Ni_map, handles.ratio_map, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Deg_2D, handles.d_Deg_2D );
else
disp('No calculated data to display')
end
disp('Done.')
guidata(hObject,handles)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('running... It may take long time to finish...')
uiwait(msgbox('Note: it may take long time to finish...'));
msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[ handles.Al_map, handles.Ni_map, handles.ratio_map ] = XY_search( handles.tot_Det_num, handles.search_Deg_2D, handles.d_Deg_2D, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_tilt2D=2;
end
if (handles.chk_cal_tilt2D ==2)
%XY_display2_solidAngle( handles.chk_print, handles.Al_map, handles.Ni_map, handles.ratio_map, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Deg_2D, handles.d_Deg_2D );
XY_display2_counts( handles.chk_print, handles.Al_map, handles.Ni_map, handles.ratio_map, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Deg_2D, handles.d_Deg_2D );
else
disp('No calculated data to display')
end
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('running... It may take long time to finish...')
uiwait(msgbox('Note: it may take long time to finish...'));
msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[ handles.Al_map_pos, handles.Ni_map_pos, handles.ratio_map_pos ] = positionXY_search( handles.tot_Det_num, handles.search_Range_2D, handles.d_Range_2D, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_shift2D=2;
end
if (handles.chk_cal_shift2D ==2)
PositionXY_display2( handles.chk_print, handles.Al_map_pos, handles.Ni_map_pos, handles.ratio_map_pos, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Range_2D, handles.d_Range_2D );
else
disp('No calculated data to display')
end
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('running... It may take long time to finish...')
uiwait(msgbox('Note: it may take long time to finish...'));
msgbox('Running...');
if (handles.chk_display == 2)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
[ handles.Al_map_pos, handles.Ni_map_pos, handles.ratio_map_pos ] = positionXY_search( handles.tot_Det_num, handles.search_Range_2D, handles.d_Range_2D, handles.sample_para, handles.holder_para, handles.holder_frame_para, handles.angle_search, handles.SpuriousX);
handles.chk_cal_shift2D=2;
end
if (handles.chk_cal_shift2D ==2)
PositionXY_display2_counts( handles.chk_print, handles.Al_map_pos, handles.Ni_map_pos, handles.ratio_map_pos, handles.tot_Det_num, handles.angle_search, handles.sample_para, handles.search_Range_2D, handles.d_Range_2D );
else
disp('No calculated data to display')
end
disp('Done.')
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)




function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
handles.search_Deg=str2double(get(hObject,'String'));
handles.search_Deg_2D=handles.search_Deg;
handles.search_Deg_2D=abs(handles.search_Deg_2D);
if (handles.search_Deg_2D>35)
    handles.search_Deg_2D=35;
    handles.search_Deg=handles.search_Deg_2D;
    disp('Warning! Our of range!')
    set(handles.edit20, 'String', num2str(handles.search_Deg_2D));
end
disp('Tilt range is set as +-');disp(handles.search_Deg_2D);
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double
handles.d_Deg=str2double(get(hObject,'String'));
handles.d_Deg_2D=handles.d_Deg;
handles.d_Deg_2D=abs(handles.d_Deg_2D);
if (handles.d_Deg_2D>handles.search_Deg_2D/3)
    handles.d_Deg_2D=handles.search_Deg_2D/3;
    handles.d_Deg=handles.d_Deg_2D;
    disp('Warning! Oversize!')
    set(handles.edit21, 'String', num2str(handles.d_Deg_2D));
end
disp('Tilt step is set as +-')
disp(handles.d_Deg_2D)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
handles.search_Range_2D=str2double(get(hObject,'String'));
handles.search_Range_2D=abs(handles.search_Range_2D);
if (handles.search_Range_2D>1)
    handles.search_Range_2D=1;
    disp('Warning! Our of range!')
    set(handles.edit22, 'String', num2str(handles.search_Range_2D));
end
disp('Shift range is set as +-')
disp(handles.search_Range_2D)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double
handles.d_Range_2D=str2double(get(hObject,'String'));
handles.d_Range_2D=abs(handles.d_Range_2D);
if (handles.d_Range_2D>handles.search_Range_2D/3)
    handles.d_Range_2D=handles.search_Range_2D/3;
    disp('Warning! Oversize!')
    set(handles.edit23, 'String', num2str(handles.d_Range_2D));
end
disp('Shift step is set as ')
disp(handles.d_Range_2D)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
chk = get(hObject,'Value');
if chk == 1 %check tiltX
    handles.chk_print=1;
    disp('Figures will be output')
else
    handles.chk_print=2;
    disp('Figures will not be output')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)

% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox11
chk = get(hObject,'Value');
if chk == 1 %check tiltX
    handles.chk_display=1;
    disp('Only display calculated results.')
else
    handles.chk_display=2;
    disp('Recalculate results when display')
end
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)


%To compiler to exe file  mcc -m SuperAngle.m



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
disp('Clear all figures')
%clear all;
h=findobj('Type','figure');
delete(h(2:length(h)));
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename_Sinput, pathname] = uigetfile({'*.xlsx'},'Update specimen information from file');
if (filename_Sinput ~=0)
%[ handles.detector_para, handles.tot_Det_num ] = Detector_input( filename );
%[ handles.angle_search, handles.detector_para] = Detector_setup( handles.detector_para, handles.dAngle );
%set(handles.staticText1, 'String', num2str(value));
%set(handles.file_det, 'String', filename);
%update_detector_info(hObject, eventdata, handles);

%demo file filename_Sinput='specimen_Ni3Al_demo_X.xlsx';
[ handles.sample_para ] = Specimen_setup( filename_Sinput, handles.t_chk );
if (handles.sample_para(12)~=1)
    set(handles.checkbox1,'Value',1);
    set(handles.edit2, 'String', num2str(handles.sample_para(17)));
end
set(handles.edit3, 'String', num2str(handles.sample_para(5)));

handles.chk_Shadow = 2;
set(handles.checkbox3,'Value',0);

set(handles.edit4, 'String', num2str(handles.sample_para(18)));
set(handles.edit5, 'String', num2str(handles.sample_para(9)));
set(handles.edit6, 'String', num2str(handles.sample_para(10)));
set(handles.edit7, 'String', num2str(handles.sample_para(21)));
set(handles.edit8, 'String', num2str(handles.sample_para(3)));
set(handles.edit9, 'String', num2str(handles.sample_para(4)));
set(handles.edit10, 'String', num2str(handles.sample_para(1)));
set(handles.edit11, 'String', num2str(handles.sample_para(2)));
set(handles.edit18, 'String', num2str(handles.sample_para(22)));
set(handles.edit19, 'String', num2str(handles.sample_para(23)));
set(handles.popupmenu1, 'Value', handles.sample_para(13)-2);
set(handles.popupmenu2, 'Value', handles.sample_para(15)-2);
set(handles.popupmenu3, 'Value', handles.sample_para(14));
set(handles.popupmenu4, 'Value', handles.sample_para(16));

if (strcmp (filename_Sinput, 'specimen_Ni3Al_demo.xlsx')) %demo mode, has Mo slot grid, setup here
    handles.holder_para(7)= 1; %grid_chk(Shadow by sample grid? 1-Yes 2-No)
    handles.holder_para(8)= 1; %type_grid(standard 1x2 slot? (1=yes others no))
    handles.holder_para(9)= 0; %open_diameter_grid(If others, input grid open diameter (in mm))
    set(handles.checkbox5,'Value',1);
    handles.exp_file='Exp_data_Ni3Al_demo_X.xlsx';
    handles.A_exp_counts=12691; %total counts
    handles.B_exp_counts=83422;
    set(handles.edit16, 'String', num2str(handles.A_exp_counts));
    set(handles.edit17, 'String', num2str(handles.B_exp_counts));
    disp('Demo Ni3Al specimen info is loaded, have fun!')
    disp ('Please find details in "A Numerical Model for Multiple Detector')
    disp ('Energy Dispersive X-ray Spectroscopy in the Transmission Electron')
    disp ('Microscope", Ultramicroscopy, 2016, in press');
    uiwait(msgbox('Demo Ni3Al specimen is loaded, have fun!'));
end
set(handles.text101, 'String', filename_Sinput);
set(handles.text67, 'String', handles.exp_file);
[ cwText ] = GetCommandWindow();
set(handles.edit24, 'String', cwText);
guidata(hObject,handles)
end
