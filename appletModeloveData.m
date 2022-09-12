function varargout = appletModeloveData(varargin)
% APPLETMODELOVEDATA M-file for appletModeloveData.fig
%      APPLETMODELOVEDATA, by itself, creates a new APPLETMODELOVEDATA or raises the existing
%      singleton*.
%
%      H = APPLETMODELOVEDATA returns the handle to a new APPLETMODELOVEDATA or the handle to
%      the existing singleton*.
%
%      APPLETMODELOVEDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPLETMODELOVEDATA.M with the given input arguments.
%
%      APPLETMODELOVEDATA('Property','Value',...) creates a new APPLETMODELOVEDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before appletModeloveData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to appletModeloveData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help appletModeloveData

% Last Modified by GUIDE v2.5 27-Jan-2013 21:18:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @appletModeloveData_OpeningFcn, ...
                   'gui_OutputFcn',  @appletModeloveData_OutputFcn, ...
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


% --- Executes just before appletModeloveData is made visible.
function appletModeloveData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to appletModeloveData (see VARARGIN)

% Choose default command line output for appletModeloveData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes appletModeloveData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = appletModeloveData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%slider dekrementu:
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
set(handles.edit5,'String',num2str(d))
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

vykonaj(A1,fm1,A2,fm2,d,handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%A1
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
set(handles.edit6,'String',num2str(A1))
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

vykonaj(A1,fm1,A2,fm2,d,handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%fm1
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.edit7,'String',num2str(fm1))
vykonaj(A1,fm1,A2,fm2,d,handles)


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%A2
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.edit8,'String',num2str(A2))
vykonaj(A1,fm1,A2,fm2,d,handles)



% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%fm2
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.edit9,'String',num2str(fm2))
vykonaj(A1,fm1,A2,fm2,d,handles)



% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function vykonaj(A1,fm1,A2,fm2,d,handles)
%vykonaj˙ sa prÌkazy öecky
t = linspace(0,1,40001);
t =t';
sig = A1*exp(-d*fm1.*t).*sin(2*pi*fm1.*t)+...
    A2*exp(-d*fm2.*t).*sin(2*pi*fm2.*t);

data(:,1) = t;
data(:,2) = sig;

save('E:\Google drive\Google Drive\Skola\doktorandura\pulzna metoda merania rezonanËnej frekvencie\analyza metodiky\zaznam_kmitania\modelData.mat','data');
N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%v˝poËet fft, 
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostrannÈ spektrum,delenie N je kÙly normovaniu asi len
Fs = N/data(end,1);%vzorkovacia frekvencia(Ëasovanie musÌ zaËÌnaù od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor prÌsluön˝ch frekvenciÌ
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
if length(M)==1
    B1 = M(1);
    f1 = f(X(I(1)));
    set(handles.edit15,'String',num2str(f1))
    set(handles.edit16,'String','Nan')
    set(handles.edit10,'String',num2str(B1))
    B2 = B1;
    f2 = f1;
else
    B1 = M(1);
    f1 = f(X(I(1)));
    B2 = M(2);
    f2 = f(X(I(2)));
    set(handles.edit15,'String',num2str(f1))
    set(handles.edit16,'String',num2str(f2))
    set(handles.edit10,'String',num2str(B1))
    set(handles.edit11,'String',num2str(B2))
    set(handles.edit12,'String',num2str(B1/B2))
end




delete(get(handles.axes1,'Children'))
hold(handles.axes1,'on')
plot(handles.axes1,t,sig);
title(handles.axes1,'Sign·l')
xlabel(handles.axes1,'Ëas / s')
ylabel(handles.axes1,'Amplit˙da / -')

delete(get(handles.axes2,'Children'))
hold(handles.axes2,'on')
plot(handles.axes2,f,spektrum)
plot(handles.axes2,f1, B1,'*','Color','red')
plot(handles.axes2,f2, B2,'*','Color','green')
title(handles.axes2,'Amplit˙dy kmitania vo frekvenËnej oblasti')
xlabel(handles.axes2,'Frekvencia / Hz')
ylabel(handles.axes2,'Amplit˙da / -')

legend(handles.axes2,'Amp. spektrum','B1', 'B2','Location','Best')
if get(handles.checkbox1,'Value')
    if fm1<fm2
        xmin = fm1-0.1*(abs(fm2-fm1));
        xmax = fm2+0.1*(abs(fm2-fm1));
        ymax = B1+B1*0.1;
        xlim(handles.axes2,[xmin xmax])
        ylim(handles.axes2,[0 ymax])
    else
        if fm1 ~= fm2
            xmin = fm2-0.05*(abs(fm2-fm1));
            xmax = fm1+0.05*(abs(fm2-fm1));
            ymax = B1+B1*0.05;
            xlim(handles.axes2,[xmin xmax])
            ylim(handles.axes2,[0 ymax])
        end
    end
else
    xmin = str2num(get(handles.edit13,'String'));
    xmax = str2num(get(handles.edit14,'String'));
    xlim(handles.axes2,[xmin xmax])
end




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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


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
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.edit5,'String');
if iscell(d)
    d = d{1};
end
d = str2num(d);
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.slider1,'Value',d)
vykonaj(A1,fm1,A2,fm2,d,handles)

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
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.edit6,'String');
A1 = str2num(A1);
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.slider2,'Value',A1)
vykonaj(A1,fm1,A2,fm2,d,handles)


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
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.edit7,'String');
if iscell(fm1)
    fm1 = fm1{1};
end
fm1 = str2num(fm1);
A2 = get(handles.slider4,'Value');
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.slider3,'Value',(fm1-fm1min)/(fm1max-fm1min))
vykonaj(A1,fm1,A2,fm2,d,handles)


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
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
A2 = get(handles.edit8,'String');
if iscell(A2)
    A2 = A2{1};
end
A2 = str2num(A2);
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;

set(handles.slider4,'Value',A2)
vykonaj(A1,fm1,A2,fm2,d,handles)


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
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
A2 = get(handles.slider4,'Value');
fm2 = get(handles.edit9,'String');
if iscell(fm2)
   fm2 = fm2{1}; 
end
fm2 = str2num(fm2);

set(handles.slider5,'Value',(fm2-fm2min)/(fm2max-fm2min))
vykonaj(A1,fm1,A2,fm2,d,handles)


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



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


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


%xmin:
function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xmin = get(handles.edit13,'String');
xmin = str2num(xmin);
xmax = get(handles.edit14,'String');
xmax = str2num(xmax);
xlim(handles.axes2,[xmin xmax])

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
xmin = get(handles.edit13,'String');
xmin = str2num(xmin);
xmax = get(handles.edit14,'String');
xmax = str2num(xmax);
xlim(handles.axes2,[xmin xmax])


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



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fm1min = str2num(get(handles.edit1,'String'));
fm1max = str2num(get(handles.edit2,'String'));
fm2min = str2num(get(handles.edit3,'String'));
fm2max = str2num(get(handles.edit4,'String'));
d = get(handles.slider1,'Value');
A1 = get(handles.slider2,'Value');
fm1 = get(handles.slider3,'Value')*(fm1max-fm1min)+fm1min;
fm2 = get(handles.slider5,'Value')*(fm2max-fm2min)+fm2min;
A2 = get(handles.slider4,'Value');

t = linspace(0,1,40001);
t =t';
%model data:
sig = A1*exp(-d*fm1.*t).*sin(2*pi*fm1.*t)+...
    A2*exp(-d*fm2.*t).*sin(2*pi*fm2.*t);

data(:,1) = t;
data(:,2) = sig;

N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%v˝poËet fft, 
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostrannÈ spektrum,delenie N je kÙly normovaniu asi len
Fs = N/data(end,1);%vzorkovacia frekvencia(Ëasovanie musÌ zaËÌnaù od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor prÌsluön˝ch frekvenciÌ
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
if length(M)==1
    B1 = M(1);
    f1 = f(X(I(1)));
    B2 = B1;
    f2 = f1;
else
    B1 = M(1);
    f1 = f(X(I(1)));
    B2 = M(2);
    f2 = f(X(I(2)));
end

figure
subplot(2,1,1)
plot(t,sig);
title('Sign·l')
xlabel('Ëas / s')
ylabel('Amplit˙da / -')

subplot(2,1,2)
plot(f,spektrum)
hold on
plot(f1, B1,'*','Color','red')
plot(f2, B2,'*','Color','green')
title('Amplit˙dy kmitania vo frekvenËnej oblasti')
xlabel('Frekvencia / Hz')
ylabel('Amplit˙da / -')
legend('Amp. spektrum','B1', 'B2','Location','Best')

if get(handles.checkbox1,'Value')
    if fm1<fm2
        xmin = fm1-0.1*(abs(fm2-fm1));
        xmax = fm2+0.1*(abs(fm2-fm1));
        ymax = B1+B1*0.1;
        xlim([xmin xmax])
        ylim([0 ymax])
    else
        if fm1 ~= fm2
            xmin = fm2-0.05*(abs(fm2-fm1));
            xmax = fm1+0.05*(abs(fm2-fm1));
            ymax = B1+B1*0.05;
            xlim([xmin xmax])
            ylim([0 ymax])
        end
    end
else
    xmin = str2num(get(handles.edit13,'String'));
    xmax = str2num(get(handles.edit14,'String'));
    xlim([xmin xmax])
end
