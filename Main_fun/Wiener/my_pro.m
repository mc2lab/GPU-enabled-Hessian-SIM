function varargout = my_pro(varargin)
%MY_PRO MATLAB code file for my_pro.fig
%      MY_PRO, by itself, creates a new MY_PRO or raises the existing
%      singleton*.
%
%      H = MY_PRO returns the handle to a new MY_PRO or the handle to
%      the existing singleton*.
%
%      MY_PRO('Property','Value',...) creates a new MY_PRO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to my_pro_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MY_PRO('CALLBACK') and MY_PRO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MY_PRO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help my_pro

% Last Modified by GUIDE v2.5 12-Dec-2017 16:30:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @my_pro_OpeningFcn, ...
                   'gui_OutputFcn',  @my_pro_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before my_pro is made visible.
function my_pro_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for my_pro
set(handles.edit1,'String',2);
set(handles.edit2,'String',488);
set(handles.edit3,'String',20);
set(handles.edit5,'String',150);
set(handles.edit6,'String',1);
set(handles.edit7,'String',65);
set(handles.edit8,'String',1.4);
set(handles.edit9,'String',4);
set(handles.edit10,'String',4);
set(handles.edit11,'String',3);
set(handles.axes1,'xTick',[]);
set(handles.axes1,'ytick',[]);
set(handles.axes1,'box','on');
set(handles.radiobutton4,'value',1);
set(handles.radiobutton6,'value',1);
set(handles.text15,'visible','off');
set(handles.text16,'visible','off');
% set(handles.axes1,'visible','off')
handles.output = hObject;

% Update handles structure
% guidata(hObject, handles);

% UIWAIT makes my_pro wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is only for Hessian-Denoise function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radiobutton2_Callback(hObject, eventdata, handles);
set(handles.radiobutton1,'visible','off');
set(handles.radiobutton3,'String','Hessian-CPU');
%set(handles.radiobutton3,'visible','off');
set(handles.radiobutton2,'String','Hessian-GPU');

% --- Outputs from this function are returned to the command line.
function varargout = my_pro_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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




function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


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


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
set(handles.radiobutton1,'value',1);
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',0);
%%
set(handles.edit1,'visible','on');
set(handles.edit2,'visible','on');
set(handles.edit3,'visible','on');
set(handles.edit5,'visible','on');
set(handles.edit6,'visible','on');
set(handles.edit7,'visible','on');
set(handles.edit8,'visible','on');
set(handles.edit9,'visible','on');
set(handles.edit10,'visible','on');
set(handles.edit11,'visible','on');
set(handles.radiobutton4,'visible','on');
set(handles.radiobutton5,'visible','on');
set(handles.radiobutton6,'visible','on');
set(handles.radiobutton7,'visible','on');
%%
set(handles.text1,'visible','on');
set(handles.text2,'visible','on');
set(handles.text3,'visible','on');
set(handles.text5,'visible','on');
set(handles.text6,'visible','on');
set(handles.text7,'visible','on');
set(handles.text8,'visible','on');
set(handles.text9,'visible','on');
set(handles.text10,'visible','on');
set(handles.text11,'visible','on');
set(handles.text12,'visible','on');
set(handles.text13,'visible','on');
set(handles.text14,'visible','on');
if get(handles.radiobutton5,'Value')
    set(handles.text16,'visible','on');
end
if get(handles.radiobutton7,'Value')
    set(handles.text15,'visible','on');
end
%%


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
set(handles.radiobutton1,'value',0);
set(handles.radiobutton2,'value',1);
set(handles.radiobutton3,'value',0);
%%
set(handles.edit1,'visible','off');
set(handles.edit2,'visible','off');
set(handles.edit3,'visible','off');
set(handles.edit5,'visible','on');
set(handles.edit6,'visible','on');
set(handles.edit7,'visible','off');
set(handles.edit8,'visible','off');
set(handles.edit9,'visible','off');
set(handles.edit10,'visible','off');
set(handles.edit11,'visible','off');
set(handles.radiobutton4,'visible','off');
set(handles.radiobutton5,'visible','off');
set(handles.radiobutton6,'visible','off');
set(handles.radiobutton7,'visible','off');
%%
set(handles.text1,'visible','off');
set(handles.text2,'visible','off');
set(handles.text3,'visible','off');
set(handles.text5,'visible','on');
set(handles.text6,'visible','on');
set(handles.text7,'visible','off');
set(handles.text8,'visible','off');
set(handles.text9,'visible','off');
set(handles.text10,'visible','on');
set(handles.text11,'visible','off');
set(handles.text12,'visible','off');
set(handles.text13,'visible','off');
set(handles.text14,'visible','off');
set(handles.text15,'visible','off');
set(handles.text16,'visible','off');
%%

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
set(handles.radiobutton2,'value',0);
set(handles.radiobutton3,'value',1);
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
% set(handles.radiobutton1,'value',0);
% set(handles.radiobutton2,'value',0);
% set(handles.radiobutton3,'value',1);
% %%
% set(handles.edit1,'visible','on');
% set(handles.edit2,'visible','on');
% set(handles.edit3,'visible','on');
% set(handles.edit5,'visible','on');
% set(handles.edit6,'visible','on');
% set(handles.edit7,'visible','on');
% set(handles.edit8,'visible','on');
% set(handles.edit9,'visible','on');
% set(handles.edit10,'visible','on');
% set(handles.edit11,'visible','on');
% set(handles.radiobutton4,'visible','on');
% set(handles.radiobutton5,'visible','on');
% set(handles.radiobutton6,'visible','on');
% set(handles.radiobutton7,'visible','on');
% %%
% set(handles.text1,'visible','on');
% set(handles.text2,'visible','on');
% set(handles.text3,'visible','on');
% set(handles.text5,'visible','on');
% set(handles.text6,'visible','on');
% set(handles.text7,'visible','on');
% set(handles.text8,'visible','on');
% set(handles.text9,'visible','on');
% set(handles.text10,'visible','on');
% set(handles.text11,'visible','on');
% set(handles.text12,'visible','on');
% set(handles.text13,'visible','on');
% set(handles.text14,'visible','on');
% if get(handles.radiobutton5,'Value')
%     set(handles.text16,'visible','on');
% end
% if get(handles.radiobutton7,'Value')
%     set(handles.text15,'visible','on');
% end
%%




% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.tif'},'Please choose the Tif format file');
if isequal(filename,0)
    warndlg('Interrupt choosing a data','warn','modal'); 
else
   disp(['User selected', fullfile(pathname, filename)])
   f_load=imreadstack(fullfile(pathname, filename));
    axes(handles.axes1);
    imshow(f_load(:,:,1),[]);
    handles.f_load=f_load;
    handles.pathname=pathname;
    handles.filename=filename;
end
guidata(hObject,handles);
% save('lastfile.mat','pathname','filename');


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if     get(handles.radiobutton1,'Value')==1 && get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton3,'Value')==0
%    weila_21_low_signal
if get(handles.radiobutton1,'Value')==0 && (get(handles.radiobutton2,'Value')==1 || get(handles.radiobutton3,'Value')==1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error checking added. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(handles, 'pathname') == 0
        helpdlg('Image is not loaded.')
        return
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section 7. Code Structure: Scripts vs. functions 
    % Make a script into a function for the MATLAB runtime to explore a reduced search space.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if handles.radiobutton2.Value == 1
        % Function call by passing parameter
        mu=str2double(get(handles.edit5,'String'));
        sigma=str2double(get(handles.edit6,'String'));
        Bregman_Hessian_Denoise(mu, sigma, handles.pathname, handles.filename)
    elseif handles.radiobutton3.Value == 1
        % Load script (original code) 
        Bregman_Hessian_Denoise_ori_script
    end 
elseif get(handles.radiobutton1,'Value')==0 && get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton3,'Value')==1
    weila_21_low_signal
elseif get(handles.radiobutton1,'Value')==0 && get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton3,'Value')==0
    warndlg('please choose a process mode','warn','modal'); 
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


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


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
set(handles.text16,'visible','off');
set(handles.radiobutton4,'value',1);
set(handles.radiobutton5,'value',0);


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
set(handles.radiobutton4,'value',0);
set(handles.radiobutton5,'value',1);
[otf_filename,otf_pathname]=uigetfile({'*.tif'},'Please choose the Tif format OTF');
if isequal(otf_filename,0)
    warndlg('Interrupt choosing a special OTF','warn','modal'); 
    set(handles.radiobutton4,'value',1);
    set(handles.radiobutton5,'value',0);
else
    disp(['User selected', fullfile(otf_pathname, otf_filename)])
    handles.otf_pathname=otf_pathname;
    handles.otf_filename=otf_filename;
    set(handles.text16,'String',[otf_pathname otf_filename]);
    set(handles.text16,'visible','on');
end
guidata(hObject,handles);


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
set(handles.text15,'visible','off');
set(handles.radiobutton6,'value',1);
set(handles.radiobutton7,'value',0);


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7
set(handles.radiobutton6,'value',0);
set(handles.radiobutton7,'value',1);
[bg_filename,bg_pathname]=uigetfile({'*.tif'},'Please choose the Tif format background');
if isequal(bg_filename,0)
    warndlg('Interrupt choosing a special background','warn','modal'); 
    set(handles.radiobutton6,'value',1);
    set(handles.radiobutton7,'value',0);
else
    disp(['User selected', fullfile(bg_pathname, bg_filename)])
    handles.bg_pathname=bg_pathname;
    handles.bg_filename=bg_filename;
    set(handles.text15,'String',[bg_pathname bg_filename]);
    set(handles.text15,'visible','on');
end
guidata(hObject,handles);
