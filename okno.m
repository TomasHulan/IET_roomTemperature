function varargout = okno(varargin)
%ako vytvori? èasovaè:
%t = 0;èas
%T = timer();%vytvorí èasovaè
%T.Period = 10;vykoná každých desa? sekúnd
%T.TimerFcn = 't+1';%zvyši premennú t o jednu sekundu
%start(T);spustí timer
%stop(T);vypne timer

% OKNO M-file for okno.fig
%      OKNO, by itself, creates a new OKNO or raises the existing
%      singleton*.
%
%      H = OKNO returns the handle to a new OKNO or the handle to
%      the existing singleton*.
%
%      OKNO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OKNO.M with the given input arguments.
%
%      OKNO('Property','Value',...) creates a new OKNO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before okno_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to okno_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help okno

% Last Modified by GUIDE v2.5 03-Jun-2015 12:30:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @okno_OpeningFcn, ...
                   'gui_OutputFcn',  @okno_OutputFcn, ...
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
end


% --- Executes just before okno is made visible.
function okno_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to okno (see VARARGIN)

% Choose default command line output for okno
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end

% UIWAIT makes okno wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = okno_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%funkcia, kde sa spracuváva nameraný signál, hladá sa dekrement a
%vykreslujú sa údaje:
function spracuj(y,t,hObject, eventdata, handles)
%y - nameraný signál
%t - trvanie nahrávania
%hObject,eventdata, handles - skopírova? z argumentu nadradenej funkcie
y = y-mean(y);%posunutie priemeru na nulu

data(:,2) = y;%uloženie dat do tabu¾ky
data(:,1) = linspace(0,t,length(data));


N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%výpoèet fft, 
if get(handles.radiobutton1,'Value')
    spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
else
    spektrum = 2*(abs(Y(1:NFFT/2+1))/N).^2;%power spectrum
end
Fs = N/data(end,1);%vzorkovacia frekvencia(èasovanie musí zaèína? od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
[A P] = ampCheck(spektrum);%hodnoty a polohy lokálnych maxím
[M m] = sort(A,'descend');%triedenie maxim zostupne
I = P(m);%pozicie maxim zostupne
index = 1;%pozicia konkrétneho maxima v postupností maxím
rezFrek = f(I(index));

%nájdenie log. dekrementu z útlmu kmitania:
[Z,X] = ampCheck(data(:,2));
cufit = fit(data(X,1),Z','exp1');%fitovamie exponencialov
koeficienty = coeffvalues(cufit);
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
logDek = -koeficienty(2)/rezFrek;

%% urèenie log. dekrementu s polšírky spekrálnej krivky:
    %lavý svah:
    peakleft = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakleft(end+1,1) = spektrum(i);
        peakleft(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i-1;
    end
    %pravý svah:
    peakright = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakright(end+1,1) = spektrum(i);
        peakright(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i+1;
    end
    if size(peakleft,1)>1 && size(peakright,1)>1
        fl = interp1(peakleft(:,1),peakleft(:,2),M(index)/2);%frekvencia na lavom svahu frek krivky v polovyške
        fr = interp1(peakright(:,1),peakright(:,2),M(index)/2);%frekvencia na pravom svahu frek krivky v polovyške
        h = fr-fl;%šírka krivky v polvýške
        logDek2 = h*pi/(sqrt(3)*rezFrek);
        set(handles.edit22,'String',num2str(logDek2));
        %iba z ¾avého svahu:
            h = 2*(rezFrek-fl);%predpokladá sa symetria
            logDekL = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit23,'String',num2str(logDekL));
        %iba z pravého svahu:
            h = 2*(fr-rezFrek);
            logDekR = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit24,'String',num2str(logDekR));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
set(handles.text11,'String',num2str(logDek));
set(handles.edit42,'String',num2str(rezFrek));
%zobrazenie výstupu:
hold(handles.axes3,'off')
sigdek = plot(handles.axes3,data(:,1),data(:,2),data(:,1),cufit(data(:,1)));%plotovanie aj s dekrementom
set(sigdek(2),'Color','green','LineWidth',1.5)%farba èiary fitovania
title(handles.axes3,'Signal y(t)')
xlabel(handles.axes3,'Time (s)')
ylabel(handles.axes3,'Voltage (V)')
set(handles.axes3,'UserData',{sigdek});

hold(handles.axes4,'off');
plot(handles.axes4,f,spektrum);
hold(handles.axes4, 'on')
plot(handles.axes4,f(I(index)),M(index),'*','Color','red');%vykreslennie hviezdièky na maxime
if size(peakleft,1)>1 && size(peakright,1)>1 %ak bol schopný nájs? polovýšku
    plot(handles.axes4,[fr fl],[M(index)/2 M(index)/2],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky v polovici svahov
else
    plot(handles.axes4,[0 0],[0 0],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky na zaèiatku sústavy
end
hold(handles.axes4,'off');
rozlisenie = f(2)-f(1);
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end
%limity frekvnèného grafu:
if get(handles.checkbox2,'Value')
    pushbutton3_Callback(hObject, eventdata, handles)
end
%vymazanie údajov o polohách uzlov:
set(handles.pushbutton24,'UserData',{});
end


%tlacidlo Start:
function pushbutton1_Callback(hObject, eventdata, handles)
%stlaèenie tlaèidla Start
set(handles.text5,'Visible','on')%zapnutie textu "runing"
pause(0.000001);%pauza aby sa zobrazil napis runing
triger = str2num(get(handles.edit3,'String'));%percentuálny nárast potrebný na vypoèet
konec = 1;%podmienka ukonèenia
%Fs = 40000;%vzorkovacia frekvencia
Fs = str2num(get(handles.edit5,'String'));
t = 0.01;%interval poèas ktorého zis?ujem èi nenastal úder
t1 = str2num(get(handles.edit4,'String'));%èas v sekundách kolko nahrávam po klepnutí
sucetold = [];
while konec
    y = wavrecord(t*Fs,Fs);
    y = y(10:end);%odseknutie zaèiatku
    if isempty(sucetold)%ak zaèinam tj. sucet old ešte neni definovany
        sucetold = sum(ampCheck(y));
        sucetnew = sum(ampCheck(y));
    else
        sucetnew = sum(ampCheck(y));%suma lokálnych maxím
    end    
    narast = (sucetnew-sucetold)/sucetold;%percentualny narast noveho suctu oproti staremu
    sucetold = sucetnew;%priradenie noveho suctu do stareho pre dalsi ciklus
    set(handles.text2,'String',num2str(narast))
    if narast>triger%triger
        y = wavrecord(t1*Fs,Fs);
        konec = 0;%ukonèenie naèúvania
    end
end
y = y-mean(y);%posunutie priemeru na nulu

data(:,2) = y;
data(:,1) = linspace(0,t1,length(data));

[Z,X] = ampCheck(data(:,2));%získanie lokálnych maxím
cufit = fit(data(X,1),Z','exp1');%fitovamie exponencialov
koeficienty = coeffvalues(cufit);
N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%výpoèet fft, 
if get(handles.radiobutton1,'Value')
    spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
else
    spektrum = 2*(abs(Y(1:NFFT/2+1))/N).^2;%power spectrum
end
Fs = N/data(end,1);%vzorkovacia frekvencia(èasovanie musí zaèína? od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
[M,I] = sort(spektrum,'descend');%triedenie amxim zostupne
index = 1;%pozicia konkrétneho maxima v postupností maxím
rezFrek = f(I(index));
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
logDek = -koeficienty(2)/rezFrek;
set(handles.text11,'String',num2str(logDek));
set(handles.edit42,'String',num2str(rezFrek));
%zobrazenie výstupu:
axes(handles.axes3)%nastavenie axes3 ako aktualnej
plot(data(:,1),data(:,2),data(:,1),cufit(data(:,1)))%plotovanie aj s dekrementom
title('Signál y(t)')
xlabel('Èas (s)')
ylabel('Napätie (V)')
axes(handles.axes4) %nastavenie axes4 ako aktualnej
plot(f,spektrum)
hold on
plot(f(I(index)),M(index),'*','Color','red')%vykreslennie hviezdièky na maxime
hold off
rozlisenie = f(2)-f(1);%najmenší dielik na frekvnènej stupnici
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end

set(handles.text5,'Visible','off')%vypnutie textu "runing"
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end


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
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end


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
end


%rozsah spektra:
function pushbutton3_Callback(hObject, eventdata, handles)
%nastavý rozsah grafu
xmin = str2num(get(handles.edit1,'String'));
xmax = str2num(get(handles.edit2,'String'));
xlim(handles.axes4,[xmin xmax])
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
end

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
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
end


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
end


%posun rez frek dopredu pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
%set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
userData = get(handles.figure1,'UserData');
index = userData{2};
index = index+1;
f = userData{3};
I = userData{1};
M = userData{5};
rezFrek = f(I(index));
set(handles.edit42,'String',num2str(rezFrek))
koeficienty = userData{4};
logDek = -koeficienty(2)/rezFrek;


    %urèenie log. dekrementu s polšírky spekrálnej krivky:
        %set(handles.pushbutton6,'UserData',{Y NFFT Fs data1 data});%uloženie fourierovho obrazu pre pásmovú zádrž
        userdata = get(handles.pushbutton6,'UserData');
        Y = userdata{1};
        NFFT = userdata{2};
        N = userdata{4};
        N = length(N);
        spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum
        
    %% urèenie log. dekrementu s polšírky spekrálnej krivky:
    %lavý svah:
    peakleft = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakleft(end+1,1) = spektrum(i);
        peakleft(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i-1;
    end
    %pravý svah:
    peakright = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakright(end+1,1) = spektrum(i);
        peakright(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i+1;
    end
    if size(peakleft,1)>1 && size(peakright,1)>1
        fl = interp1(peakleft(:,1),peakleft(:,2),M(index)/2);%frekvencia na lavom svahu frek krivky v polovyške
        fr = interp1(peakright(:,1),peakright(:,2),M(index)/2);%frekvencia na pravom svahu frek krivky v polovyške
        h = fr-fl;%šírka krivky v polvýške
        logDek2 = h*pi/(sqrt(3)*rezFrek);
        set(handles.edit22,'String',num2str(logDek2));
        %iba z ¾avého svahu:
            h = 2*(rezFrek-fl);%predpokladá sa symetria
            logDekL = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit23,'String',num2str(logDekL));
        %iba z pravého svahu:
            h = 2*(fr-rezFrek);
            logDekR = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit24,'String',num2str(logDekR));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
hold(handles.axes4,'off');
plot(handles.axes4,f,spektrum);
hold(handles.axes4, 'on')
plot(handles.axes4,f(I(index)),M(index),'*','Color','red');%vykreslennie hviezdièky na maxime
if size(peakleft,1)>1 && size(peakright,1)>1 %ak bol schopný nájs? polovýšku
    plot(handles.axes4,[fr fl],[M(index)/2 M(index)/2],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky v polovici svahov
else
    plot(handles.axes4,[0 0],[0 0],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky na zaèiatku sústavy
end
hold(handles.axes4,'off');
%limity frekvnèného grafu:
if get(handles.checkbox2,'Value')
    pushbutton3_Callback(hObject, eventdata, handles)
end
rozlisenie = f(2)-f(1);
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end
end


%dozadu:
function pushbutton5_Callback(hObject, eventdata, handles)
%set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
userData = get(handles.figure1,'UserData');
index = userData{2};
if index>1%ak som neni na zaèiatku
index = index-1;
f = userData{3};
I = userData{1};
M = userData{5};
rezFrek = f(I(index));
set(handles.edit42,'String',num2str(rezFrek))
koeficienty = userData{4};
logDek = -koeficienty(2)/rezFrek;
    %urèenie log. dekrementu s polšírky spekrálnej krivky:
        %set(handles.pushbutton6,'UserData',{Y NFFT Fs data1 data});%uloženie fourierovho obrazu pre pásmovú zádrž
        userdata = get(handles.pushbutton6,'UserData');
        Y = userdata{1};
        NFFT = userdata{2};
        N = userdata{4};
        N = length(N);
        spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum
        
    %% urèenie log. dekrementu s polšírky spekrálnej krivky:
    %lavý svah:
    peakleft = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakleft(end+1,1) = spektrum(i);
        peakleft(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i-1;
    end
    %pravý svah:
    peakright = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        peakright(end+1,1) = spektrum(i);
        peakright(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
        i = i+1;
    end
    if size(peakleft,1)>1 && size(peakright,1)>1
        fl = interp1(peakleft(:,1),peakleft(:,2),M(index)/2);%frekvencia na lavom svahu frek krivky v polovyške
        fr = interp1(peakright(:,1),peakright(:,2),M(index)/2);%frekvencia na pravom svahu frek krivky v polovyške
        h = fr-fl;%šírka krivky v polvýške
        logDek2 = h*pi/(sqrt(3)*rezFrek);
        set(handles.edit22,'String',num2str(logDek2));
        %iba z ¾avého svahu:
            h = 2*(rezFrek-fl);%predpokladá sa symetria
            logDekL = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit23,'String',num2str(logDekL));
        %iba z pravého svahu:
            h = 2*(fr-rezFrek);
            logDekR = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit24,'String',num2str(logDekR));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
hold(handles.axes4,'off');
plot(handles.axes4,f,spektrum);
hold(handles.axes4, 'on')
plot(handles.axes4,f(I(index)),M(index),'*','Color','red');%vykreslennie hviezdièky na maxime
if size(peakleft,1)>1 && size(peakright,1)>1 %ak bol schopný nájs? polovýšku
    plot(handles.axes4,[fr fl],[M(index)/2 M(index)/2],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky v polovici svahov
else
    plot(handles.axes4,[0 0],[0 0],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky na zaèiatku sústavy
end
hold(handles.axes4,'off');
rozlisenie = f(2)-f(1);
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end
%limity frekvnèného grafu:
if get(handles.checkbox2,'Value')
    pushbutton3_Callback(hObject, eventdata, handles)
end
end
end



function edit5_Callback(hObject, eventdata, handles)
%nic
end

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
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
end

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
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
end

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
end


%pasmova zadrž:
function pushbutton6_Callback(hObject, eventdata, handles)

end


%stlaèenie tlaèidla "Pulz":
function pushbutton7_Callback(hObject, eventdata, handles)
pause(0.0000001);%pauza aby sa zobrazil napis runing
set(handles.text5,'Visible','on')%zapnutie textu "runing"
amplituda = str2num(get(handles.edit8,'String'));
casImpulzu = str2num(get(handles.edit17,'String'));
frekPulzu = str2num(get(handles.edit20,'String'));%frekvencia, ktorá je posielaná poèas pulzu
%Fs = 40000;%vzorkovacia frekvencia
Fs = str2num(get(handles.edit5,'String'));
t1 = str2num(get(handles.edit4,'String'));%èas v sekundách kolko nahrávam po klepnutí
% if amplituda >1
%     set(handles.edit8,'String','1')
%     amplituda = 1;
% end
if amplituda<0
    set(handles.edit8,'String','0')
    amplituda = 0;
end
t = 0:0.0001:casImpulzu;%
pulz = amplituda*sin(2*pi*frekPulzu.*t);%èíslo je frekvencia pulzu
wavplay(pulz,round(length(t)/casImpulzu));%impulz
pause(str2num(get(handles.edit9,'String')))%oneskorenie po vyslaní signálu
y = wavrecord(t1*Fs,Fs);%nahra? zvuk
spracuj(y,t1,hObject, eventdata, handles)%spracovnie dát a zobrazenie výstupov
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
%vypnutie textu "runing"
set(handles.text5,'Visible','off')
end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
end

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
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
end

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
end



function edit10_Callback(hObject, eventdata, handles)
end

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
end


function edit11_Callback(hObject, eventdata, handles)
end

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
end

% Tlaèidlo na orezávanie:
function pushbutton8_Callback(hObject, eventdata, handles)
userData = get(handles.pushbutton6,'UserData');%naèitanie dat
data(:,1) = userData{4};
data(:,2) = userData{5};
od = str2num(get(handles.edit10,'String'));%odreza? od (v sekundách)
do = str2num(get(handles.edit11,'String'));%odreza? do (v sekundách)
if do > data(end,1)%ak je horna hranica mimo nastavý sa najvyšší èas
   do = data(end,1); 
end
%najdenie inexov na odrezanie:
interval = (data(2,1)-data(1,1))/2;
iod = find(data(:,1)>=od-interval,1,'first');
ido = find(data(:,1)>=do-interval,1,'first');
data1(:,1) = data(iod:ido,1)-data(iod,1);%naèitanie x-ovych dat a posunutie na nulu
data1(:,2) = data(iod:ido,2);
data1(:,2) = data1(:,2)-mean(data1(:,2));%posunutie priemeru na nulu
data = data1;
%spracovanie a zobrazenie:
spracuj(data(:,2),data(end,1),hObject, eventdata, handles)
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

%oreza? z¾ava:
function pushbutton10_Callback(hObject, eventdata, handles)
userData = get(handles.pushbutton6,'UserData');%naèitanie dat
data(:,1) = userData{4};
data(:,2) = userData{5};
krok = str2num(get(handles.edit12,'String'));
%orezanie z¾ava o "krok" udajov:
data1(:,1) = data(krok+1:end,1)-data(krok+1,1);
data1(:,2) = data(krok+1:end,2);
data1(:,2) = data1(:,2)-mean(data1(:,2));%posunutie priemeru na nulu
%spracovanie a zobrazenie dát
spracuj(data1(:,2),data1(end,1),hObject, eventdata, handles)
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
end

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
end



%pásmová priepus?:
function pushbutton11_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
cas = userdata{4};
od = str2num(get(handles.edit15,'String'));
do = str2num(get(handles.edit16,'String'));
nh = round(do*NFFT/Fs);%najvyšší index
nd = round(od*NFFT/Fs);%najnižší index
maska = [zeros(nd,1); ones(nh-nd,1); zeros(NFFT-2*nh,1);...
    ones(nh-nd,1); zeros(nd,1)];
Y = Y.*maska;%aplikácia pásmovej priepuste
data = ifft(Y,NFFT);
data = real(data(1:length(cas)));%vybranie realneho priebehu po filtrovani
spracuj(data,cas(end),hObject, eventdata, handles)
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
end

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
end


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
end

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
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
end

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
end

function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
end

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
end

function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
end

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
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
end

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
end


% priepus? okolo rezonanènej frekvencie
function pushbutton15_Callback(hObject, eventdata, handles)
userData = get(handles.figure1,'UserData');
%set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
userData1 = get(handles.pushbutton6,'UserData');
I = userData{1};%zostupné poradie pozícií najväèších amplitúd
index = userData{2};%pozícia aktuálne selektovanej frekvencie v premennej I
f = userData{3};%x-ová os (frekvencie)
Y = userData1{1};%komplexné spektrum
cas = userData1{4};
NFFT = userData1{2};
Fs = userData1{3};
N = length(cas);%poèet nameraných dát
sirka = str2num(get(handles.edit21,'String'));%šírka pasma prepustených frekvencií okolo rezonanènej
vyber1 = f<f(I(index))+sirka/2;%podmienka pre hornú hranicu
vyber2 = f>f(I(index))-sirka/2;%podmienaka pre dolnú hranicu
vyber = vyber1.*vyber2;%logický prienik podmienok pre hornu a dolnu hranicu
rvyber = [vyber vyber(end-2:-1:1)];%rozšírený výber pre komplexné spektrum

Y = Y.*rvyber';%aplikácia masky(filtra) na dvojstranné komplexné spektrum

data = ifft(Y,NFFT);
data = real(data(1:N));%vybranie realneho priebehu po filtrovani
%spracovanie a zobrazenie vystupu:
spracuj(data,cas(end),hObject, eventdata, handles)
end


function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double
end

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
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
end


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
axes(handles.axes3) %nastavenie axes4 ako aktualnej
if strcmp(get(gca,'YGrid'),'off')
    set(gca,'YGrid','on')
else
    set(gca,'YGrid','off')
end
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
end

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
end

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
%set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
userdata = get(handles.pushbutton6,'UserData');
data(:,1) = userdata{4};
data(:,2) = userdata{5};
%uloži? do špecifického miesta:
subor = get(handles.edit29,'String');
cesta = get(handles.edit30,'String');
save([cesta subor '.mat'],'data');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double
end

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
end


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double
end

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
end


function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double
end

% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double
end

% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double
end

% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%Nahra?
function pushbutton19_Callback(hObject, eventdata, handles)
%naèítanie dát zo súboru:
    subor = get(handles.edit28,'String');
    cesta = get(handles.edit30,'String');
    y = load([cesta subor '.mat']);
    y = y.data;%vybra? data zo štruktúry
%%%%%%%%%%%%%%%%%%%%%%%%%
%spracovanie a vykreslenie:
spracuj(y(:,2),y(end,1)-y(1,1),hObject, eventdata, handles)
end


function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
end

% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double
end

% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double
end

% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double
end

% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double
end

% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
%set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));
%f1 musí by? menšie ako f2 ak nie tak sa vymenia:
if f1>f2
    pom = f1;
    f1 = f2;
    f2 = pom;
    pom = A1;
    A1 = A2;
    A2 = pom;
end
set(handles.edit31,'String',num2str(f1));
set(handles.edit32,'String',num2str(f2));

%nájdenie maxím signálu:
[Z,X] = localMax(data(:,2),3);
%nájdenie minimá z maxím (uzly signálu):
Z = smooth(Z);
[Zmin,Xmin] = localMin(Z,10);
% [M,I] = sort(Z,'descend');%triedenie maxim zostupne
indNearMax = round((X(Xmin(2))+X(Xmin(1)))/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
[Zmax,Xmax] = localMax(data(indNearMax:end,2),3);%haldám maximá už len zo zvyšku
indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
T = Zmax(1);%suèet amplitúd v modeli, T = B1+B2; aj ked nie presne, lebo je závislé od dekremntu -> nelinearita
hold(handles.axes3,'on')
plot(handles.axes3,[data(X(Xmin(1)),1) data(X(Xmin(2)),1)], [Zmin(1) Zmin(2)],'diamond','Color','red')
plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');



%set(handles.axes4,'UserData',{gspektrum greddot ggreendot});%uloženie handle na grafy
userdata = get(handles.axes4,'UserData');
greddot = userdata{2};
delete(greddot)
hold(handles.axes4,'on');
greddot = plot(handles.axes4,[f1 f2],[A1 A2],'*','Color','red');
userdata{2} = greddot;
set(handles.axes4,'UserData',userdata);%uloženie handle na grafy
hold(handles.axes4,'off')

matica = [1 -A1*f1/(A2*f2);1 1];
vektor = [0;T];
B = inv(matica)*vektor;
set(handles.edit33,'String',num2str(B(1)));
set(handles.edit34,'String',num2str(B(2)));
set(handles.edit36,'String',num2str(B(1)/B(2)));


%nájdenie bodu napravo od maxima, ktorý má nulovú výchylku:
i = 0;
konec = 1;
while konec
    i = i+1;
    if data(indMax-i,2) <= 0
        indZero = indMax-i;%index bodu, ktorý  prešiel nulou
        konec = 0;
    end
end
%cas, ktorému prislúcha nulová výchylka (aby fitovanie zaèínalo s bodom s
%nulovou výchylkou):
zeroPoint = interp1([data(indZero,2) data(indZero+1,2)], [data(indZero,1) data(indZero+1,1)], 0);
cas = [zeroPoint; data(indZero:end,1)];
cas = cas-zeroPoint;%aby zaèínal od nuly(kôly fitovaniu)
vychylky = [0; data(indZero:end,2)];
[Z,X] = localMax(data(:,2),1);
maxAmp = max(Z);%maximálna amplitúda dosiahnutá v nahratom signále
%nájdenie bodu, ktorý dosihol hornú hranicu rozsahu zaznamenaného signálu
%(pri nahrávaní mikrofónom je max amplitúda 1), alebo ma maximálnu výchylku:
for i = indZero:-1:1
    if (data(i,2))>=maxAmp
        firstPoint = i;
        break
    end
end
minuscas = data(firstPoint:indZero,1)-zeroPoint;%casy pred nulou
minuleVychylky = data(firstPoint:indZero,2);%vychylky prisluchajuce casom pred t=0s.
cas = [minuscas; cas];
vychylky = [minuleVychylky; vychylky];

%fitovanie ja zadam amplitudy::    
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,...
               'Upper',1,...
               'Startpoint',0.362);
%odhaduje sa parameter "a" èo je logartimický dekrement
g = fittype('B1*exp(-a*f1.*x).*sin(2*pi*f1.*x)+B2*exp(-a*f2.*x).*sin(2*pi*f2.*x)','problem',{'B1' 'f1' 'B2' 'f2'},'options',s);
[c1,gof1] = fit(cas,vychylky,g,'problem',{B(1) f1 B(2) f2});
coef = coeffvalues(c1);%vypèítaný(odhadnutý parameter) dekrement útlmu
set(handles.edit41,'String',num2str(coef(1)))
set(handles.edit50,'String',num2str(coef(1)))
set(handles.text49,'String',num2str(gof1.rsquare));
y = c1(cas);
hold(handles.axes3,'on')
[upy iupy] = localMax(y,1);
[downy idowny] = localMin(y,1);
plot(handles.axes3,cas(iupy)+zeroPoint,upy,'Color','yellow');
plot(handles.axes3,cas(idowny)+zeroPoint,downy,'Color','yellow');
hold(handles.axes3,'off')


%fitovanie fitujú sa aj amilitúdy:
s = fitoptions('Method','NonlinearLeastSquares');
s.Lower = [0 0 0];
s.Upper = [1 100 100];
s.Startpoint = [0.687 0.511 0.911];
s.Robust = 'Off';
s.Algorithm = 'Trust-Region';
s.DiffMinChange = 1e-8;
s.DiffMaxChange = 0.1;
s.MaxFunEvals = 600;
s.MaxIter = 400;
s.TolFun = 1e-6;
s.TolX = 1e-6;
g = fittype('b*exp(-a*f1.*x).*sin(2*pi*f1.*x)+c*exp(-a*f2.*x).*sin(2*pi*f2.*x)','problem',{'f1' 'f2'},'options',s);

[c2,gof2] = fit(cas,vychylky,g,'problem',{f1 f2});
coef = coeffvalues(c2);%vypèítaný(odhadnutý parameter) dekrement útlmu
d = coef(1);%dekrement
set(handles.edit38,'String',num2str(coef(2)))
set(handles.edit39,'String',num2str(coef(3)))
set(handles.edit40,'String',num2str(coef(2)/coef(3)))
set(handles.edit35,'String',num2str(d))
set(handles.edit51,'String',num2str(d))
set(handles.text42,'String',num2str(gof2.rsquare));
y = c2(cas);
hold(handles.axes3,'on')
[upy iupy] = localMax(y,1);
[downy idowny] = localMin(y,1);
plot(handles.axes3,cas(iupy)+zeroPoint,upy,'Color','red');
plot(handles.axes3,cas(idowny)+zeroPoint,downy,'Color','red');
hold(handles.axes3,'off')

grafy = get(handles.axes3,'Children');
for i = 1:length(grafy)
    xudaje = get(grafy(i),'XData');
    set(grafy(i),'Xdata',xudaje-zeroPoint)
end
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double
end

% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
end

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double
end

% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double
end

% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double
end

% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double
end

% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double
end

% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double
end

% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double
end

% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%new fig
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
%po kliknutí na toto tlaèidlo sa objavý nová figúra s dvoma grafmi
chax3 = get(handles.axes3,'Children');
chax4 = get(handles.axes4,'Children');
nchax3 = length(chax3);
nchax4 = length(chax4);

figure
%
subplot(2,1,1)
xpopis = (get(handles.axes3,'XLabel'));
xlabel(get(xpopis,'String'))
ypopis = (get(handles.axes3,'YLabel'));
ylabel(get(ypopis,'String'))
titul = (get(handles.axes3,'Title'));
title(get(titul,'String'))
hold on
for i = nchax3:-1:1
    xdata = get(chax3(i),'XData');
    ydata = get(chax3(i),'YData');
    linestyle = get(chax3(i),'LineStyle');
    linewidth = get(chax3(i),'LineWidth');
    marker = get(chax3(i),'Marker');
    color = get(chax3(i),'Color');  
    plot(xdata,ydata,'LineStyle',linestyle,'LineWidth',linewidth...
        ,'Marker',marker,'Color',color);
end
hold off
%
subplot(2,1,2)
xpopis = (get(handles.axes4,'XLabel'));
xlabel(get(xpopis,'String'))
ypopis = (get(handles.axes4,'YLabel'));
ylabel(get(ypopis,'String'))
titul = (get(handles.axes4,'Title'));
title(get(titul,'String'))
hold on
for i = nchax4:-1:1
    xdata = get(chax4(i),'XData');
    ydata = get(chax4(i),'YData');
    linestyle = get(chax4(i),'LineStyle');
    linewidth = get(chax4(i),'LineWidth');
    marker = get(chax4(i),'Marker');
    color = get(chax4(i),'Color');  
    plot(xdata,ydata,'LineStyle',linestyle,'LineWidth',linewidth...
        ,'Marker',marker,'Color',color);
end
xlim(get(handles.axes4,'XLim'))
hold off
%
end


%rucne hladanie uzlov <
function pushbutton21_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));

%%%%
userdata = get(handles.axes4,'UserData');
greddot = userdata{2};
delete(greddot)
hold(handles.axes4,'on');
greddot = plot(handles.axes4,[f1 f2],[A1 A2],'*','Color','red');
userdata{2} = greddot;
set(handles.axes4,'UserData',userdata);%uloženie handle na grafy
hold(handles.axes4,'off')
%%%%
posun = str2num(get(handles.edit43,'String'));
perioda = 1/abs(f1-f2);%perioda rázov
br = perioda*Fs;%nameraných bodov na jeden ráz
hold(handles.axes3,'on')
userdata = get(handles.pushbutton24,'UserData');
if isempty(userdata)==0
    p1 = userdata{1};
    p2 = userdata{2};
    p3 = userdata{3};
    u1 = userdata{4};
    u2 = userdata{5};
end
if isempty(userdata)
    p1=plot(handles.axes3,[0.1 0.1], [-1 1],'Color','black');
    p2=plot(handles.axes3,[0.1+perioda 0.1+perioda], [-1 1],'Color','black');
    p3=plot(handles.axes3,[0.1+2*perioda 0.1+2*perioda], [-1 1],'Color','black');
    u1 = 0.1*Fs;
    u2 = u1+br;
    set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
    indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
    [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
    indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
    hold(handles.axes3,'on')
    userdata = get(handles.axes3,'UserData');
    if length(userdata)>1
        delete(userdata{2})%vymazanie bodky v strede rázu
    end
    greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
    set(handles.axes3,'UserData',{userdata{1} greddot})
else  
    if u1-posun >0
        u1 = u1-posun;
        u2 = u2-posun;
        if ishandle(p1)
            delete(p1);delete(p2);delete(p3);
        end
        p1=plot(handles.axes3,[u1/Fs u1/Fs], [-1 1],'Color','black');
        p2=plot(handles.axes3,[u1/Fs+perioda u1/Fs+perioda], [-1 1],'Color','black');
        p3=plot(handles.axes3,[u1/Fs+2*perioda u1/Fs+2*perioda], [-1 1],'Color','black');
        set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
            indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
        [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
        indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
        hold(handles.axes3,'on')
        userdata = get(handles.axes3,'UserData');
        if length(userdata)>1
            delete(userdata{2})%vymazanie bodky v strede rázu
        end
        greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
        set(handles.axes3,'UserData',{userdata{1} greddot})
    end
end
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));
userdata = get(handles.pushbutton24,'UserData');
if isempty(userdata)==0
    p1 = userdata{1};
    p2 = userdata{2};
    p3 = userdata{3};
    u1 = userdata{4};
    u2 = userdata{5};
end
posun = str2num(get(handles.edit43,'String'));
perioda = 1/abs(f1-f2);%perioda rázov
br = perioda*Fs;%nameraných bodov na jeden ráz
hold(handles.axes3,'on')
if isempty(userdata)
    p1=plot(handles.axes3,[0.1 0.1], [-1 1],'Color','black');
    p2=plot(handles.axes3,[0.1+perioda 0.1+perioda], [-1 1],'Color','black');
    p3=plot(handles.axes3,[0.1+2*perioda 0.1+2*perioda], [-1 1],'Color','black');
    u1 = 0.1*Fs;
    u2 = u1+br;
    set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
    indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
    [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
    indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
    hold(handles.axes3,'on')
    userdata = get(handles.axes3,'UserData');
    if length(userdata)>1
        delete(userdata{2})%vymazanie bodky v strede rázu
    end
    greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
    set(handles.axes3,'UserData',{userdata{1} greddot})
else  
    if u1-posun >0
        u1 = u1+posun;
        u2 = u2+posun;
        if ishandle(p1)
            delete(p1);delete(p2);delete(p3);
        end
        p1=plot(handles.axes3,[u1/Fs u1/Fs], [-1 1],'Color','black');
        p2=plot(handles.axes3,[u1/Fs+perioda u1/Fs+perioda], [-1 1],'Color','black');
        p3=plot(handles.axes3,[u1/Fs+2*perioda u1/Fs+2*perioda], [-1 1],'Color','black');
        set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
        indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
        [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
        indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
        hold(handles.axes3,'on')
        userdata = get(handles.axes3,'UserData');
        if length(userdata)>1
            delete(userdata{2})%vymazanie bodky v strede rázu
        end
        greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
        set(handles.axes3,'UserData',{userdata{1} greddot})
    end
end
end

%tlaèidlo "Ruène"
function pushbutton24_Callback(hObject, eventdata, handles)
%set(handles.pushbutton6,'UserData',{Y NFFT Fs data(:,1) data(:,2)});%uloženie fourierovho obrazu pre pásmovú zádrž
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));
%f1 musí by? menšie ako f2 ak nie tak sa vymenia:
if f1>f2
    pom = f1;
    f1 = f2;
    f2 = pom;
    pom = A1;
    A1 = A2;
    A2 = pom;
end
set(handles.edit31,'String',num2str(f1));
set(handles.edit32,'String',num2str(f2));

%set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
userdata = get(handles.pushbutton24,'UserData');
u1 = userdata{4};
u2 = userdata{5};
indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
[Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
T = Zmax(1);%suèet amplitúd v modeli, T = B1+B2; aj ked nie presne, lebo je závislé od dekremntu -> nelinearita
hold(handles.axes3,'on')
plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');



%set(handles.axes4,'UserData',{gspektrum greddot ggreendot});%uloženie handle na grafy
userdata = get(handles.axes4,'UserData');
greddot = userdata{2};
delete(greddot)
hold(handles.axes4,'on');
greddot = plot(handles.axes4,[f1 f2],[A1 A2],'*','Color','red');
userdata{2} = greddot;
set(handles.axes4,'UserData',userdata);%uloženie handle na grafy
hold(handles.axes4,'off')

matica = [1 -A1*f1/(A2*f2);1 1];
vektor = [0;T];
B = inv(matica)*vektor;
set(handles.edit33,'String',num2str(B(1)));
set(handles.edit34,'String',num2str(B(2)));
set(handles.edit36,'String',num2str(B(1)/B(2)));


%nájdenie bodu napravo od maxima, ktorý má nulovú výchylku:
i = 0;
konec = 1;
while konec
    i = i+1;
    if data(indMax-i,2) <= 0
        indZero = indMax-i;%index bodu, ktorý má prešiel nulou
        konec = 0;
    end
end
%cas, ktorému prislúcha nulová výchylka (aby fitovanie zaèínalo s bodom s
%nulovou výchylkou):
zeroPoint = interp1([data(indZero,2) data(indZero+1,2)], [data(indZero,1) data(indZero+1,1)], 0);
cas = [zeroPoint; data(indZero:end,1)];
cas = cas-zeroPoint;%aby zaèínal od nuly(kôly fitovaniu)
vychylky = [0; data(indZero:end,2)];
[Z,X] = localMax(data(:,2),1);
maxAmp = max(Z);%maximálna amplitúda dosiahnutá v nahratom signále
%nájdenie bodu, ktorý dosihol hornú hranicu rozsahu zaznamenaného signálu
%(pri nahrávaní mikrofónom je max amplitúda 1), alebo ma maximálnu výchylku:
for i = indZero:-1:1
    if (data(i,2))>=maxAmp
        firstPoint = i;
        break
    end
end
minuscas = data(firstPoint:indZero,1)-zeroPoint;%casy pred nulou
minuleVychylky = data(firstPoint:indZero,2);%vychylky prisluchajuce casom pred t=0s.
cas = [minuscas; cas];
vychylky = [minuleVychylky; vychylky];

%fitovanie ja zadam amplitudy::    
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,...
               'Upper',1,...
               'Startpoint',0.362);
%odhaduje sa parameter "a" èo je logartimický dekrement
g = fittype('B1*exp(-a*f1.*x).*sin(2*pi*f1.*x)+B2*exp(-a*f2.*x).*sin(2*pi*f2.*x)','problem',{'B1' 'f1' 'B2' 'f2'},'options',s);
[c1,gof1] = fit(cas,vychylky,g,'problem',{B(1) f1 B(2) f2});
coef = coeffvalues(c1);%vypèítaný(odhadnutý parameter) dekrement útlmu
set(handles.edit41,'String',num2str(coef(1)))
set(handles.text49,'String',num2str(gof1.rsquare));
y = c1(cas);
hold(handles.axes3,'on')
[upy iupy] = localMax(y,1);
[downy idowny] = localMin(y,1);
plot(handles.axes3,cas(iupy)+zeroPoint,upy,'Color','yellow');
plot(handles.axes3,cas(idowny)+zeroPoint,downy,'Color','yellow');
hold(handles.axes3,'off')


%fitovanie fitujú sa aj amilitúdy:
s = fitoptions('Method','NonlinearLeastSquares');
s.Lower = [0 0 0];
s.Upper = [1 100 100];
s.Startpoint = [0.687 0.511 0.911];
s.Robust = 'Off';
s.Algorithm = 'Trust-Region';
s.DiffMinChange = 1e-8;
s.DiffMaxChange = 0.1;
s.MaxFunEvals = 600;
s.MaxIter = 400;
s.TolFun = 1e-6;
s.TolX = 1e-6;
g = fittype('b*exp(-a*f1.*x).*sin(2*pi*f1.*x)+c*exp(-a*f2.*x).*sin(2*pi*f2.*x)','problem',{'f1' 'f2'},'options',s);

[c2,gof2] = fit(cas,vychylky,g,'problem',{f1 f2});
coef = coeffvalues(c2);%vypèítaný(odhadnutý parameter) dekrement útlmu
d = coef(1);%dekrement
set(handles.edit38,'String',num2str(coef(2)))
set(handles.edit39,'String',num2str(coef(3)))
set(handles.edit40,'String',num2str(coef(2)/coef(3)))
set(handles.edit35,'String',num2str(d));
set(handles.text42,'String',num2str(gof2.rsquare));
y = c2(cas);
hold(handles.axes3,'on')
[upy iupy] = localMax(y,1);
[downy idowny] = localMin(y,1);
plot(handles.axes3,cas(iupy)+zeroPoint,upy,'Color','red');
plot(handles.axes3,cas(idowny)+zeroPoint,downy,'Color','red');
hold(handles.axes3,'off')
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double
end

% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function pushbutton15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


% --- Executes during object creation, after setting all properties.
function pushbutton8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end



function radiobutton1_ButtonDownFcn(hObject, eventdata, handles)

end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
data(:,1) = userdata{4};
data(:,2) = userdata{5};

[Z,X] = ampCheck(data(:,2));
cufit = fit(data(X,1),Z','exp1');%fitovamie exponencialov
koeficienty = coeffvalues(cufit);
N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%výpoèet fft, 
if get(handles.radiobutton1,'Value')
    spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
else
    spektrum = 2*(abs(Y(1:NFFT/2+1))/N).^2;%power spectrum
end
Fs = N/data(end,1);%vzorkovacia frekvencia(èasovanie musí zaèína? od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[M,I] = sort(spektrum,'descend');%triedenie maxim zostupne
index = 1;%pozicia konkrétneho maxima v postupností maxím
rezFrek = f(I(index));
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
logDek = -koeficienty(2)/rezFrek;

%urèenie log. dekrementu s polšírky spekrálnej krivky:
    %lavý svah:
    peakleft = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        i = i-1;
        peakleft(end+1,1) = spektrum(i);
        peakleft(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
    end
    %pravý svah:
    peakright = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        i = i+1;
        peakright(end+1,1) = spektrum(i);
        peakright(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
    end
    if size(peakleft,1)>1 && size(peakright,1)>1
        fl = interp1(peakleft(:,1),peakleft(:,2),M(index)/2);%frekvencia na lavom svahu frek krivky v polovyške
        fr = interp1(peakright(:,1),peakright(:,2),M(index)/2);%frekvencia na pravom svahu frek krivky v polovyške
        h = fr-fl;%šírka krivky v polvýške
        logDek2 = h*pi/(sqrt(3)*rezFrek);
        set(handles.edit22,'String',num2str(logDek2));
        %iba z ¾avého svahu:
            h = 2*(rezFrek-fl);%predpokladá sa symetria
            logDekL = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit23,'String',num2str(logDekL));
        %iba z pravého svahu:
            h = 2*(fr-rezFrek);
            logDekR = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit24,'String',num2str(logDekR));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.text11,'String',num2str(logDek));
set(handles.edit42,'String',num2str(rezFrek));
%zobrazenie výstupu:
hold(handles.axes3,'off')
sigdek = plot(handles.axes3,data(:,1),data(:,2),data(:,1),cufit(data(:,1)));%plotovanie aj s dekrementom
title(handles.axes3,'Signal y(t)')
xlabel(handles.axes3,'Time (s)')
ylabel(handles.axes3,'Voltage (V)')
set(handles.axes3,'UserData',{sigdek});

hold(handles.axes4,'off');
gspektrum = plot(handles.axes4,f,spektrum);
hold(handles.axes4, 'on')
greddot = plot(handles.axes4,f(I(index)),M(index),'*','Color','red');%vykreslennie hviezdièky na maxime
if size(peakleft,1)>1 && size(peakright,1)>1 %ak bol schopný nájs? polovýšku
    ggreendot = plot(handles.axes4,[fr fl],[M(index)/2 M(index)/2],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky v polovici svahov
else
    ggreendot = plot(handles.axes4,[0 0],[0 0],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky na zaèiatku sústavy
end
hold(handles.axes4,'off');
set(handles.axes4,'UserData',{gspektrum greddot ggreendot});%uloženie handle na grafy
%xlim([0 Fs*t1])%zobrazi? frekvencie iba do 10000 Hz
rozlisenie = f(2)-f(1);
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end
set(handles.text5,'Visible','off')%vypnutie textu "runing"
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
%limity frekvnèného grafu:
if get(handles.checkbox2,'Value')
    pushbutton3_Callback(hObject, eventdata, handles)
end
%vymazanie údajov o polohách uzlov:
set(handles.pushbutton24,'UserData',{});
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
data(:,1) = userdata{4};
data(:,2) = userdata{5};

[Z,X] = ampCheck(data(:,2));
cufit = fit(data(X,1),Z','exp1');%fitovamie exponencialov
koeficienty = coeffvalues(cufit);
N = length(data);
NFFT = 2^nextpow2(N);
Y = fft(data(:,2),NFFT);%výpoèet fft, 
if get(handles.radiobutton1,'Value')
    spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
else
    spektrum = 2*(abs(Y(1:NFFT/2+1))/N).^2;%power spectrum
end
Fs = N/data(end,1);%vzorkovacia frekvencia(èasovanie musí zaèína? od nuly)
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[M,I] = sort(spektrum,'descend');%triedenie maxim zostupne
index = 1;%pozicia konkrétneho maxima v postupností maxím
rezFrek = f(I(index));
set(handles.figure1,'UserData',{I index f koeficienty M});%uloženie triedených frekvencií I, index - konkrétne maximum
logDek = -koeficienty(2)/rezFrek;

%urèenie log. dekrementu s polšírky spekrálnej krivky:
    %lavý svah:
    peakleft = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        i = i-1;
        peakleft(end+1,1) = spektrum(i);
        peakleft(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
    end
    %pravý svah:
    peakright = [];
    konec = 1;
    i = I(index);%index rez frekvencie
    while konec
        i = i+1;
        peakright(end+1,1) = spektrum(i);
        peakright(end,2) = f(i);
        %ukonèi? ak prejdem polvýškou:
        if spektrum(i)<M(index)/2
            konec = 0;
        end
    end
    if size(peakleft,1)>1 && size(peakright,1)>1
        fl = interp1(peakleft(:,1),peakleft(:,2),M(index)/2);%frekvencia na lavom svahu frek krivky v polovyške
        fr = interp1(peakright(:,1),peakright(:,2),M(index)/2);%frekvencia na pravom svahu frek krivky v polovyške
        h = fr-fl;%šírka krivky v polvýške
        logDek2 = h*pi/(sqrt(3)*rezFrek);
        set(handles.edit22,'String',num2str(logDek2));
        %iba z ¾avého svahu:
            h = 2*(rezFrek-fl);%predpokladá sa symetria
            logDekL = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit23,'String',num2str(logDekL));
        %iba z pravého svahu:
            h = 2*(fr-rezFrek);
            logDekR = h*pi/(sqrt(3)*rezFrek);
            set(handles.edit24,'String',num2str(logDekR));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.text11,'String',num2str(logDek));
set(handles.edit42,'String',num2str(rezFrek));
%zobrazenie výstupu:
hold(handles.axes3,'off')
sigdek = plot(handles.axes3,data(:,1),data(:,2),data(:,1),cufit(data(:,1)));%plotovanie aj s dekrementom
title(handles.axes3,'Signal y(t)')
xlabel(handles.axes3,'Time (s)')
ylabel(handles.axes3,'Voltage (V)')
set(handles.axes3,'UserData',{sigdek});

hold(handles.axes4,'off');
gspektrum = plot(handles.axes4,f,spektrum);
hold(handles.axes4, 'on')
greddot = plot(handles.axes4,f(I(index)),M(index),'*','Color','red');%vykreslennie hviezdièky na maxime
if size(peakleft,1)>1 && size(peakright,1)>1 %ak bol schopný nájs? polovýšku
    ggreendot = plot(handles.axes4,[fr fl],[M(index)/2 M(index)/2],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky v polovici svahov
else
    ggreendot = plot(handles.axes4,[0 0],[0 0],'*','Color',[0.582 0.3867 0.3867]);%vykreslenie hviezdièky na zaèiatku sústavy
end
hold(handles.axes4,'off');
set(handles.axes4,'UserData',{gspektrum greddot ggreendot});%uloženie handle na grafy
%xlim([0 Fs*t1])%zobrazi? frekvencie iba do 10000 Hz
rozlisenie = f(2)-f(1);
if get(handles.radiobutton1,'Value')
    title(handles.axes4,['Single-Sided Amplitude Spetrum of y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'|Y(f)|')
else
    title(handles.axes4,['Power Spectrum of signalu y(t), resolution: ' num2str(rozlisenie) ' Hz'])
    xlabel(handles.axes4,'Frequency (Hz)')
    ylabel(handles.axes4,'Power')
end
set(handles.text5,'Visible','off')%vypnutie textu "runing"
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
%limity frekvnèného grafu:
if get(handles.checkbox2,'Value')
    pushbutton3_Callback(hObject, eventdata, handles)
end
%vymazanie údajov o polohách uzlov:
set(handles.pushbutton24,'UserData',{});
end


%urèenie dekrementu z amplitúdy rázov
function pushbutton25_Callback(hObject, eventdata, handles)

E1 = str2num(get(handles.edit46,'String'));
E2 = str2num(get(handles.edit47,'String'));
f1 = str2num(get(handles.edit44,'String'));
f2 = str2num(get(handles.edit45,'String'));
%výpoèet dekrementu (podla Nakutis, 2010):
d = (f2-f1)/mean([f1 f2])*log(E1/E2);
set(handles.edit48,'String',num2str(d));
end


function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double
end

% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double
end

% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double
end

% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double
end

% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double
end

% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));

set(handles.edit44,'String',num2str(f1));
set(handles.edit45,'String',num2str(f2));

%%%%
userdata = get(handles.axes4,'UserData');
greddot = userdata{2};
delete(greddot)
hold(handles.axes4,'on');
greddot = plot(handles.axes4,[f1 f2],[A1 A2],'*','Color','red');
userdata{2} = greddot;
set(handles.axes4,'UserData',userdata);%uloženie handle na grafy
hold(handles.axes4,'off')
%%%%
posun = str2num(get(handles.edit49,'String'));
perioda = 1/abs(f1-f2);%perioda rázov
br = round(perioda*Fs);%nameraných bodov na jeden ráz
hold(handles.axes3,'on')
userdata = get(handles.pushbutton24,'UserData');
if isempty(userdata)==0
    p1 = userdata{1};
    p2 = userdata{2};
    p3 = userdata{3};
    u1 = userdata{4};
    u2 = userdata{5};
end
if isempty(userdata)
    p1=plot(handles.axes3,[0.1 0.1], [-1 1],'Color','black');
    p2=plot(handles.axes3,[0.1+perioda 0.1+perioda], [-1 1],'Color','black');
    p3=plot(handles.axes3,[0.1+2*perioda 0.1+2*perioda], [-1 1],'Color','black');
    u1 = 0.1*Fs;
    u2 = u1+br;
    set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
    indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
    [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku    
    indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
    set(handles.edit46,'Strin',num2str(data(indMax,2)));
    %%%%%%%%%%%%%%%%%%%
        [Zmax2,Xmax2] = localMax(data(indNearMax+br:end,2),1);%hladam max druheho razu
        indMax2 = Xmax2(1)+indNearMax+br-1;%èasový index maxima druhého rázu
        set(handles.edit47,'Strin',num2str(data(indMax2,2)));
    %%%%%%%%%%%%%%%%%%%
    hold(handles.axes3,'on')
    userdata = get(handles.axes3,'UserData');
    if length(userdata)>1
        delete(userdata{2})%vymazanie bodky v strede rázu
    end
    greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
    set(handles.axes3,'UserData',{userdata{1} greddot})
else  
    if u1-posun >0
        u1 = u1-posun;
        u2 = u2-posun;
        if ishandle(p1)
            delete(p1);delete(p2);delete(p3);
        end
        p1=plot(handles.axes3,[u1/Fs u1/Fs], [-1 1],'Color','black');
        p2=plot(handles.axes3,[u1/Fs+perioda u1/Fs+perioda], [-1 1],'Color','black');
        p3=plot(handles.axes3,[u1/Fs+2*perioda u1/Fs+2*perioda], [-1 1],'Color','black');
        set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
            indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
        [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
        indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
        set(handles.edit46,'Strin',num2str(data(indMax,2)));
        %%%%%%%%%%%%%%%%%%%
            [Zmax2,Xmax2] = localMax(data(indNearMax+br:end,2),1);%hladam max druheho razu
            indMax2 = Xmax2(1)+indNearMax+br-1;%èasový index maxima druhého rázu
            set(handles.edit47,'Strin',num2str(data(indMax2,2)));
        %%%%%%%%%%%%%%%%%%%
        hold(handles.axes3,'on')
        userdata = get(handles.axes3,'UserData');
        if length(userdata)>1
            delete(userdata{2})%vymazanie bodky v strede rázu
        end
        greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
        set(handles.axes3,'UserData',{userdata{1} greddot})
    end
end

end

% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
data(:,1) = userdata{4};
data(:,2) = userdata{5};

%nájdenie frekvencií:
N = length(data);
spektrum = 2*abs(Y(1:NFFT/2+1)/N);%jednostranné spektrum,delenie N je kôly normovaniu asi len
f = Fs/2*linspace(0,1,NFFT/2+1);%priestor príslušných frekvencií
[Z,X] = ampCheck(spektrum);
[M,I] = sort(Z,'descend');%triedenie maxim zostupne
A1 = M(1);
A2 = M(2);
f1 = f(X(I(1)));
f2 = f(X(I(2)));
userdata = get(handles.pushbutton24,'UserData');
if isempty(userdata)==0
    p1 = userdata{1};
    p2 = userdata{2};
    p3 = userdata{3};
    u1 = userdata{4};
    u2 = userdata{5};
end
posun = str2num(get(handles.edit49,'String'));
perioda = 1/abs(f1-f2);%perioda rázov
br = round(perioda*Fs);%nameraných bodov na jeden ráz
hold(handles.axes3,'on')
if isempty(userdata)
    p1=plot(handles.axes3,[0.1 0.1], [-1 1],'Color','black');
    p2=plot(handles.axes3,[0.1+perioda 0.1+perioda], [-1 1],'Color','black');
    p3=plot(handles.axes3,[0.1+2*perioda 0.1+2*perioda], [-1 1],'Color','black');
    u1 = 0.1*Fs;
    u2 = u1+br;
    set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
    indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
    [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
    indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
    set(handles.edit46,'Strin',num2str(data(indMax,2)));
    %%%%%%%%%%%%%%%%%%%
        [Zmax2,Xmax2] = localMax(data(indNearMax+br:end,2),1);%hladam max druheho razu
        indMax2 = Xmax2(1)+indNearMax+br-1;%èasový index maxima druhého rázu
        set(handles.edit47,'Strin',num2str(data(indMax2,2)));
    %%%%%%%%%%%%%%%%%%%
    hold(handles.axes3,'on')
    userdata = get(handles.axes3,'UserData');
    if length(userdata)>1
        delete(userdata{2})%vymazanie bodky v strede rázu
    end
    greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
    set(handles.axes3,'UserData',{userdata{1} greddot})
else  
    if u1-posun >0
        u1 = u1+posun;
        u2 = u2+posun;
        if ishandle(p1)
            delete(p1);delete(p2);delete(p3);
        end
        p1=plot(handles.axes3,[u1/Fs u1/Fs], [-1 1],'Color','black');
        p2=plot(handles.axes3,[u1/Fs+perioda u1/Fs+perioda], [-1 1],'Color','black');
        p3=plot(handles.axes3,[u1/Fs+2*perioda u1/Fs+2*perioda], [-1 1],'Color','black');
        set(handles.pushbutton24,'UserData',{p1 p2 p3 u1 u2});
        indNearMax = round((u1+u2)/2);%index blízko maxima - prvé maximum napravo je to, v stred rázu
        [Zmax,Xmax] = localMax(data(indNearMax:end,2),1);%haldám maximá už len zo zvyšku
        indMax = Xmax(1)+indNearMax-1;%èasový index maxima najvyššieho rázu
        set(handles.edit46,'Strin',num2str(data(indMax,2)));
        %%%%%%%%%%%%%%%%%%%
            [Zmax2,Xmax2] = localMax(data(indNearMax+br:end,2),1);%hladam max druheho razu
            indMax2 = Xmax2(1)+indNearMax+br-1;%èasový index maxima druhého rázu
            set(handles.edit47,'Strin',num2str(data(indMax2,2)));
        %%%%%%%%%%%%%%%%%%%
        hold(handles.axes3,'on')
        userdata = get(handles.axes3,'UserData');
        if length(userdata)>1
            delete(userdata{2})%vymazanie bodky v strede rázu
        end
        greddot=plot(handles.axes3,data(indMax,1),data(indMax,2),'*','Color','red');
        set(handles.axes3,'UserData',{userdata{1} greddot})
    end
end
end


function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double
end

% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
end


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display('fdsaf')
end


% zobrazene--- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
lim = get(handles.axes3,'XLim');
od = lim(1);%oreza? od
do = lim(2);%oreza? do
userData = get(handles.pushbutton6,'UserData');%naèitanie dat
data(:,1) = userData{4};
data(:,2) = userData{5};
if do > data(end,1)%ak je horna hranica mimo nastavý sa najvyšší èas
   do = data(end,1); 
end
%najdenie inexov na odrezanie:
interval = (data(2,1)-data(1,1))/2;
iod = find(data(:,1)>=od-interval,1,'first');
ido = find(data(:,1)>=do-interval,1,'first');
data1(:,1) = data(iod:ido,1)-data(iod,1);%naèitanie x-ovych dat a posunutie na nulu
data1(:,2) = data(iod:ido,2);
data1(:,2) = data1(:,2)-mean(data1(:,2));%posunutie priemeru na nulu
%spracovanie a zobrazenie:
spracuj(data1(:,2),data1(end,1),hObject, eventdata, handles)
%aplikácia pásmovej priepuste:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
end

function y = nacuvac(handles)
%funkcia na "naèúvanie" kde sa vykonáva sluèka
set(handles.pushbutton29,'String','Stop')%zmena tlaèidla na "stop"
set(handles.pushbutton29,'BackgroundColor','red')
th = str2num(get(handles.edit3,'String'));%treshold
Fs = str2num(get(handles.edit5,'String'));%vzorkovacia frekvencia
t1 = str2num(get(handles.edit4,'String'));%èas v sekundách kolko nahrávam po klepnutí

yes = 0;%podmienka ukonèenia cyklu
while ~yes
    %y = wavrecord(t1*Fs,Fs);%nahrat zvuk
    recObj = audiorecorder(Fs, 8, 1);%vytvori objekt
    recordblocking(recObj, t1);%nahra zvuk
    y = getaudiodata(recObj);%priradi do premennej
    y = y-mean(y);%posunutie priemeru na nulu
    [Z,X] = ampCheck(y);
    n = length(find(Z>th));%celkový poèet prekmitov cez prahovú hodnotu th
    if n>3%ak je poèet prekmitov viacej ako tri
        yes = 1;
    end
    stop = get(handles.pushbutton29,'UserData');
    if ~stop %ukonèenie pomocou tlaèidla
        return
    end
    set(handles.text2,'String',num2str(max(Z)))
    pause(0.0001)
end
end

%Trigger - pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
stop = get(handles.pushbutton29,'UserData');
if isempty(stop)|| ~stop
    stop = 1;
else
    stop = 0;
end
set(handles.pushbutton29,'UserData',stop);%nastavenie príznaku na "aktívny"
if ~stop%ukonèi? tlaèidlom
    set(handles.pushbutton29,'String','Start')
    set(handles.pushbutton29,'BackgroundColor',[0.941 0.941 0.941])
    return
end
y = nacuvac(handles);
stop = get(handles.pushbutton29,'UserData');
if ~stop%ukonèi? tlaèidlom
    return
end
set(handles.pushbutton29,'UserData',0);%nastavenie príznaku na "neaktívny"
set(handles.pushbutton29,'String','Start')
set(handles.pushbutton29,'BackgroundColor',[0.941 0.941 0.941])

t1 = str2double(get(handles.edit4,'String'));%èas v sekundách kolko nahrávam po klepnutí
spracuj(y,t1,hObject, eventdata, handles)%spracovanie a vystup
%aplikácia pásmovej priepuste okolo rez frek:
if get(handles.checkbox1,'Value')
    pushbutton15_Callback(hObject, eventdata, handles)
end
%aplikácia pásmovej priepuste:
if get(handles.checkbox3,'Value')
    pushbutton11_Callback(hObject, eventdata, handles)
end
end



%save - uloženie údajov do tabulky
function pushbutton31_Callback(hObject, eventdata, handles)
rezF = str2double(get(handles.edit42,'String'));%rezonanèná frekvencia
rezF = sprintf('%.1f', rezF);
rezF(rezF=='.') = ',';
if get(handles.radiobutton3,'Value')==1
    dekr = get(handles.text11,'String');
    dekr(dekr=='.')=',';
end
if get(handles.radiobutton4,'Value')==1
    dekr = get(handles.edit22,'String');
    dekr(dekr=='.')=',';
end
if get(handles.radiobutton5,'Value')==1
    dekr = get(handles.edit23,'String');
    dekr(dekr=='.')=',';
end
if get(handles.radiobutton6,'Value')==1
    dekr = get(handles.edit24,'String');
    dekr(dekr=='.')=',';
end
if get(handles.radiobutton7,'Value')==1
    dekr = get(handles.edit50,'String');
    dekr(dekr=='.')=',';
end
if get(handles.radiobutton8,'Value')==1
    dekr = get(handles.edit51,'String');
    dekr(dekr=='.')=',';
end
data = get(handles.uitable2,'Data');
data{end+1,1}=rezF;
data{end,2} = dekr;
set(handles.uitable2,'Data',data);
end


%clear - vymaže tabu¾ku
function pushbutton32_Callback(hObject, eventdata, handles)
data = {};
set(handles.uitable2,'Data',data);
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double
end

% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double
end

% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%Del - pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
%vymaže posledný riadok tabu¾ky
data = get(handles.uitable2,'Data');
data=data(1:end-1,:);
set(handles.uitable2,'Data',data);
end

% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
end


%pásmová zádrž pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
userdata = get(handles.pushbutton6,'UserData');
Y = userdata{1};
NFFT = userdata{2};
Fs = userdata{3};
cas = userdata{4};
N = length(cas);
od = str2num(get(handles.edit6,'String'));
do = str2num(get(handles.edit7,'String'));
nh = round(do*NFFT/Fs);%najvyšší index
nd = round(od*NFFT/Fs);%najnižší index
Y(nd+1:nh) = 0;%aplikácia pasmovej zádrže
Y(end-nh:end-nd) = 0;
data = ifft(Y,NFFT);
data = real(data(1:length(cas)));%vybranie realneho priebehu po filtrovani
spracuj(data,cas(end),hObject, eventdata, handles)
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
end