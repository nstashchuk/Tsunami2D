function varargout = Tsunami2D(varargin)
% TSUNAMI2D MATLAB code for Tsunami2D.fig
%      TSUNAMI2D, by itself, creates a new TSUNAMI2D or raises the existing
%      singleton*.
%
%      H = TSUNAMI2D returns the handle to a new TSUNAMI2D or the handle to
%      the existing singleton*.
%
%      TSUNAMI2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TSUNAMI2D.M with the given input arguments.
%
%      TSUNAMI2D('Property','Value',...) creates a new TSUNAMI2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tsunami2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tsunami2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tsunami2D

% Last Modified by GUIDE v2.5 03-Feb-2020 15:44:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tsunami2D_OpeningFcn, ...
                   'gui_OutputFcn',  @Tsunami2D_OutputFcn, ...
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


% --- Executes just before Tsunami2D is made visible.
function Tsunami2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tsunami2D (see VARARGIN)

% Choose default command line output for Tsunami2D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tsunami2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tsunami2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function shelfdepth_value_Callback(hObject, eventdata, handles)
% hObject    handle to shelfdepth_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shelfdepth_value as text
%        str2double(get(hObject,'String')) returns contents of shelfdepth_value as a double
ShelfDepth = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function shelfdepth_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shelfdepth_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beginningofshelf_value_Callback(hObject, eventdata, handles)
% hObject    handle to beginningofshelf_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beginningofshelf_value as text
%        str2double(get(hObject,'String')) returns contents of beginningofshelf_value as a double
BeginningOfShelf = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function beginningofshelf_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beginningofshelf_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endofshelf_value_Callback(hObject, eventdata, handles)
% hObject    handle to endofshelf_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endofshelf_value as text
%        str2double(get(hObject,'String')) returns contents of endofshelf_value as a double
EndOfShelf = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function endofshelf_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endofshelf_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_topography.
function plot_topography_Callback(hObject, eventdata, handles)
% hObject    handle to plot_topography (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Lx 
clear global H Jtop
global H Jtop
Lx=2000000;                                 % Length of model domain (m)
H1=500;                                     % Depth in the open part (m)
ShelfDepth = str2double(get(handles.shelfdepth_value,'String'));
BeginningOfShelf = str2double(get(handles.beginningofshelf_value,'String'));
EndOfShelf = str2double(get(handles.endofshelf_value,'String'));
DX=str2double(get(handles.dx_value,'String'));
dx=DX*1000;
x=[0:dx:Lx];                                % Grid vector 
Jtop=length(x);                               % Number of grid points
H2=ShelfDepth;
Lx1=BeginningOfShelf*1000;
Lx2=EndOfShelf*1000;
            
             %  Block2: Definition of the bottom 
for j=1:Jtop                              % Loop for bottom definition
    H(j)=H1;                            % Deep part
    if x(j)>Lx1 & x(j)<=Lx2             % Beginning of the continental slope 
        H(j)=H1-(H1-H2)*(x(j)-Lx1)/(Lx2-Lx1);   % Linear topography from Lx1 to Lx2
    end                                 % End of the continental slope 
    if x(j)>Lx2                         % Start of the shelf
        H(j)=H2;                        % Shelf
    end                                 % End of the shelf
end             % End of the loop for bottom definition 

axes(handles.axes_forplot);

plot(x/1000,H,'LineWidth',6,'Color',[0.5 0.5 0.5]) %plot bottom 
set(gca,'YDir','reverse');
set(gca,'XTick',[0:100:Lx/1000]);
ylim([0 H1]);                          % Limits for y-axis  
grid on  
xlabel('Distance (km)');                % Create x-label
ylabel('Depth (m)');                    % Create y-label



function amplitude_value_Callback(hObject, eventdata, handles)
% hObject    handle to amplitude_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of amplitude_value as text
%        str2double(get(hObject,'String')) returns contents of amplitude_value as a double
Amplitude = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function amplitude_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to amplitude_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epicenter_value_Callback(hObject, eventdata, handles)
% hObject    handle to epicenter_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epicenter_value as text
%        str2double(get(hObject,'String')) returns contents of epicenter_value as a double
Epicenter = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function epicenter_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epicenter_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_value_Callback(hObject, eventdata, handles)
% hObject    handle to width_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_value as text
%        str2double(get(hObject,'String')) returns contents of width_value as a double
Width = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function width_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotElevation.
function PlotElevation_Callback(hObject, eventdata, handles)
% hObject    handle to PlotElevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global zz Jel
global Lx 
global zz Jel
DX=str2double(get(handles.dx_value,'String'));
dx=DX*1000;
x=[0:dx:Lx];                                % Grid vector 
Jel=length(x);
Amplitude = str2double(get(handles.amplitude_value,'String'));
Epicenter = str2double(get(handles.epicenter_value,'String'));
Width = str2double(get(handles.width_value,'String'));
Steepness = str2double(get(handles.steepness_value,'String'));
  Lc=Epicenter*1000;
  La=Width*1000;
  A=Amplitude;
  M=Steepness*2;
             Li=Lc-La/2;
             nstart=round(Li/dx)+1;
             zz=zeros(1,Jel);
for j=nstart:Jel                              % Loop for initial conditions 
    if (j-nstart+1)*dx < La                    % Initial free surface displacement
        zz(1,j)=A*(sin(pi*(j-nstart+1)*dx/La))^M;  % in the centre of the earthquake 
    end
end
% cla(handles.axes_forplot) 
axes(handles.axes_elevation);
 
plot(x/1000,zz(1,:),'LineWidth',2)  % plot elevation
set(gca,'XTick',[0:100:Lx/1000],'YTick',[-A:A/5:A]); 
% set(gca,'XMinorTick','on')
grid on                                  % Activation of a grid
ylim([-A A]);                            % Limits for y-axis    
xlabel('Distance (km)');                 % Create x-label
ylabel('Depth (m)');          



function steepness_value_Callback(hObject, eventdata, handles)
% hObject    handle to steepness_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of steepness_value as text
%        str2double(get(hObject,'String')) returns contents of steepness_value as a double
Steepness = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function steepness_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to steepness_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx_value_Callback(hObject, eventdata, handles)
% hObject    handle to dx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx_value as text
%        str2double(get(hObject,'String')) returns contents of dx_value as a double
dx = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function dx_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_value_Callback(hObject, eventdata, handles)
% hObject    handle to dt_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_value as text
%        str2double(get(hObject,'String')) returns contents of dt_value as a double
dt = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function dt_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_value_Callback(hObject, eventdata, handles)
% hObject    handle to time_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_value as text
%        str2double(get(hObject,'String')) returns contents of time_value as a double
Time = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function time_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_model.
function run_model_Callback(hObject, eventdata, handles)
% hObject    handle to run_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Lx
global H Jtop
global zz Jel
global z Nt x
DT=str2double(get(handles.dt_value,'String'));
DX=str2double(get(handles.dx_value,'String'));
Time=str2double(get(handles.time_value,'String'));
Latitude=str2double(get(handles.latitude_value,'String'));
NonLinearity=str2double(get(handles.nonlinearity_value,'String'));
fi=Latitude;
dt=DT;
Nt=round(Time*3600/dt);
dx=DX*1000;
x=[0:dx:Lx];                                % Grid vector 
Jx=length(x);
g=9.8;                                      % acceleration due to gravity                

Omega=0.0000729;                            % Angular velocity of the Earth rotation (1/rad)
K=NonLinearity;
f=2*Omega*sind(fi);
if Jx ~=Jtop
cla(handles.axes_forplot);    
axes(handles.axes_forplot);
text(0.1,0.1,'Recalculate topography','FontSize',14,'FontWeight','bold','Color','r')
elseif Jx ~= Jel
cla(handles.axes_elevation) 
axes(handles.axes_elevation);
text(0.1,0.1,'Recalculate elevation','FontSize',14,'FontWeight','bold','Color','r')
else    
v=zeros(Nt,Jx); u=v; z=v;
z(1,1:Jx)=zz(1,1:Jx);
n=1;       % first time step
    for j=2:Jx-1                      % Loop for u,v,z at first time dtep
       u(n+1,j)=u(n,j)+f*dt*v(n,j)-g*dt*(z(n,j+1)-z(n,j-1))/2/dx;
        v(n+1,j)=v(n,j)-f*dt*u(n,j);
       z(n+1,j)=z(n,j)-dt*(H(j+1)*u(n,j+1)-H(j-1)*u(n,j-1))/2/dx;
    end                          % end of the loop
%%%%%%%%%%%%%%%%% Boundary conditions at first temporal step    
    u(n+1,1)=0;                          % In the epicentre j=1
    v(n+1,1)=v(n+1,2);                   % In the epicentre j=1
    z(n+1,1)=z(n+1,2);                   % In the epicentre j=1 
     u(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx
     v(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx
     z(n+1,Jx)=0;                         % Right boundary on the shelf j=Jx    

H1=H(1);
H2=H(end);
%%%%%%%%%%%%%%%%%%% Block7: The rest of time steps (2:Nt) CTCS
 for n=2:Nt-1    
    for j=2:Jx-1
       u(n+1,j)=u(n-1,j)+f*2*dt*v(n,j)-g*dt*(z(n,j+1)-z(n,j-1))/dx-K*dt*u(n,j)*(u(n,j+1)-u(n,j-1))/dx;
       v(n+1,j)=v(n-1,j)-f*2*dt*u(n,j)-K*dt*u(n,j)*(v(n,j+1)-v(n,j-1))/dx;
       z(n+1,j)=z(n-1,j)-dt*(H(j+1)*u(n,j+1)-H(j-1)*u(n,j-1))/dx -K*dt*(z(n,j+1)*u(n,j+1)-z(n,j-1)*u(n,j-1))/dx ;
    end
 %%%%%%%%%%%%%%%%% Boundary Conditions for u,v,z

     u(n+1,Jx)=u(n,Jx)-(dt/dx)*sqrt(g*H2)*(u(n,Jx)-u(n,Jx-1));                       
     v(n+1,Jx)=v(n,Jx)-(dt/dx)*sqrt(g*H2)*(v(n,Jx)-v(n,Jx-1));                       
     z(n+1,Jx)=z(n,Jx)-(dt/dx)*sqrt(g*H2)*(z(n,Jx)-z(n,Jx-1)); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     u(n+1,1)=u(n,1)+(dt/dx)*sqrt(g*H1)*(u(n,2)-u(n,1));                       
     v(n+1,1)=v(n,1)+(dt/dx)*sqrt(g*H1)*(v(n,2)-v(n,1));                       
     z(n+1,1)=z(n,1)+(dt/dx)*sqrt(g*H1)*(z(n,2)-z(n,1));                       
 end     
  cla(handles.axes_elevation) 
axes(handles.axes_elevation);
text(0.1,0.1,'Movie is ready','FontSize',14,'FontWeight','bold')    
end    
   


% --- Executes on button press in run_movie_elevation.
function run_movie_elevation_Callback(hObject, eventdata, handles)
% hObject    handle to run_movie_elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Lx
global z Nt x 
Amplitude = str2double(get(handles.amplitude_value,'String'));
A=Amplitude;
cla(handles.axes_elevation) 
axes(handles.axes_elevation);  
Amax=max(max(abs(z)));
for i=1:5:Nt
plot(x/1000,z(i,:),'LineWidth',2)  % plot elevation
set(gca,'XTick',[0:100:Lx/1000],'YTick',[-Amax:Amax/5:Amax]); 
set(gca,'XMinorTick','on')
grid on                                  % Activation of a grid
ylim([-Amax Amax]);                            % Limits for y-axis    
xlabel('Distance (km)');                 % Create x-label
ylabel('Depth (m)');                     % Create y-label
drawnow
end



function latitude_value_Callback(hObject, eventdata, handles)
% hObject    handle to latitude_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latitude_value as text
%        str2double(get(hObject,'String')) returns contents of latitude_value as a double
Latitude = str2double(get(hObject,'String'));



% --- Executes during object creation, after setting all properties.
function latitude_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latitude_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonlinearity_value_Callback(hObject, eventdata, handles)
% hObject    handle to nonlinearity_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonlinearity_value as text
%        str2double(get(hObject,'String')) returns contents of nonlinearity_value as a double
NonLinearity = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function nonlinearity_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonlinearity_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function plot_e_Callback(hObject, eventdata, handles)
% hObject    handle to plot_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes_elevationfig = figure;
% Copy the axes and size it to the figure

axes_elevationcopy = copyobj(handles.axes_elevation,axes_elevationfig);
set(axes_elevationcopy,'Units','Normalized',...
              'Position',[.12,.12,.78,.78]);
guidata(hObject,handles);

% --------------------------------------------------------------------
function plot_ele_Callback(hObject, eventdata, handles)
% hObject    handle to plot_ele (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plot_t_Callback(hObject, eventdata, handles)
% hObject    handle to plot_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes_topofig = figure;
% Copy the axes and size it to the figure

axes_topocopy = copyobj(handles.axes_forplot,axes_topofig);
set(axes_topocopy,'Units','Normalized',...
              'Position',[.12,.12,.78,.78]);

% Save handles to new fig and axes in case
%  we want to do anything else to them
% handles.axes_spiral1fig = axes_spiral1fig;
% handles.axes_spiral1copy = axes_spiral1copy;
guidata(hObject,handles);


% --------------------------------------------------------------------
function plot_topo_Callback(hObject, eventdata, handles)
% hObject    handle to plot_topo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
