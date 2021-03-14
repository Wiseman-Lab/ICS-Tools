function varargout = sv(varargin)
% SV M-file for sv.fig
%      SV, by itself, creates a new SV or raises the existing
%      singleton*.
%
%      H = SV returns the handle to a new SV or the handle to
%      the existing singleton*.
%
%      SV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SV.M with the given input arguments.
%
%      SV('Property','Value',...) creates a new SV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sv_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sv

% Last Modified by GUIDE v2.5 11-May-2010 16:08:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sv_OpeningFcn, ...
                   'gui_OutputFcn',  @sv_OutputFcn, ...
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


% --- Executes just before sv is made visible.
function sv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sv (see VARARGIN)

% Choose default command line output for sv
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sv wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.in=length(varargin)-1; %0 if one input, 1 if 2 inputs
handles.sc = 0;
handles.colorcm = 'pink';
handles.stack1 = varargin{1};
handles.totim = size(handles.stack1,3);
if handles.in==0
 handles.mode='c'; %contour plot
 handles.ch=0;

elseif handles.in==1
 if strcmp(varargin{2},'c')  %contour plot
     handles.mode='c';
 elseif strcmp(varargin{2},'s')  %surface plot
     handles.mode='s';
     handles.maxz=max(handles.stack1(:));
     handles.minz=min(handles.stack1(:));
 end
 handles.ch=0;
 
elseif handles.in==2
 if strcmp(varargin{2},'c')  %contour plot
     handles.mode='c';
 elseif strcmp(varargin{2},'s')  %surface plot
     handles.mode='s';
     handles.maxz=max(handles.stack1(:));
     handles.minz=min(handles.stack1(:));
 end
 %3rd argument controls the speed of the display
handles.playspeed=varargin{3};
 handles.ch=0;
elseif handles.in==3
    if strcmp(varargin{3},'c')  %contour plot
     handles.mode='c';
    elseif strcmp(varargin{3},'s')  %surface plot
     handles.mode='s';
     handles.maxz=max(handles.stack1(:));
     handles.minz=min(handles.stack1(:));
    elseif strcmp(varargin{3},'2c')  %2color  
    handles.mode='2c';
    handles.stack2 = varargin{2};
    for it=1:handles.totim
    c3 = zeros(size(handles.stack2,1),size(handles.stack2,2));
     im1 = imadjust(normalisation(handles.stack1(:,:,it)));
     im2 = imadjust(normalisation(handles.stack2(:,:,it)));
     handles.out(:,:,1,it) = im1;
     handles.out(:,:,2,it) = im2;
     handles.out(:,:,3,it) = c3;
    end
    handles.out1=squeeze(handles.out(:,:,:,1));
    end
    handles.ch=1;
handles.stack2 = varargin{2};% get image data from input
handles.image2 = handles.stack2(:,:,1);  
handles.maxz2=max(handles.stack2(:));
handles.minz2=min(handles.stack2(:));
handles.playspeed=varargin{4};
end
handles.image1 = handles.stack1(:,:,1);

set(handles.imtot, 'String', handles.totim);    % diplay total number of images
% set(handles.imslider,'Min',1)
% set(handles.imslider,'Max',handles.totim)
set(gcf,'Name','Stack Viewer')
showim(handles)
set(handles.imslider,'SliderStep',[handles.playspeed*0.01/handles.totim handles.playspeed*0.1/handles.totim])
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = sv_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function imslider_Callback(hObject, eventdata, handles)
% hObject    handle to imslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
image_number = get(hObject,'Value');
image_number = round(image_number*handles.totim);
if image_number == 0 
	image_number = 1;
end
set(handles.imnb, 'String', image_number);
handles.image1 = handles.stack1(:,:,image_number);
if handles.ch==1
handles.image2 = handles.stack2(:,:,image_number);
end
showim(handles)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function imslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in sccolor.
function sccolor_Callback(hObject, eventdata, handles)
% hObject    handle to sccolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sccolor
showim(handles)
handles.sc = get(hObject,'Value');
guidata(hObject, handles);

function showim(handles)
if handles.sc == 0 
    
    if handles.in==0
    imagesc(handles.image1)
    axis on, colormap(handles.colorcm)
    
   % axis on, axis image, colormap(handles.colorcm)
    elseif handles.in==1 || handles.in==2
      if strcmp(handles.mode,'c')  
          imagesc(handles.image1)
          axis on, colormap(handles.colorcm)
      elseif strcmp(handles.mode,'s')
          surf(double(handles.image1))
          set(gca,'view',[95.5           64])%[-0.5000   50.0000]);%[ -78.5000   78.0000]) 
           axis on, colormap(handles.colorcm)
           zlim([handles.minz handles.maxz])
      end
    elseif handles.in==3
      if strcmp(handles.mode,'c') 
          subplot(1,2,1)
          imagesc(handles.image1)
          subplot(1,2,2)
          imagesc(handles.image2)
          axis on, colormap(handles.colorcm)
      elseif strcmp(handles.mode,'s')
          subplot(1,2,1)
          surf(double(handles.image1))
          set(gca,'view',[-0.5000   50.0000]);%[ -78.5000   78.0000]) 
          zlim([handles.minz handles.maxz])
          axis on, colormap(handles.colorcm)
          subplot(1,2,2)
          surf(double(handles.image2))
          zlim([handles.minz2 handles.maxz2])
          set(gca,'view',[-0.5000   50.0000]);%[ -78.5000   78.0000]) 
          axis on, colormap(handles.colorcm) 
      elseif strcmp(handles.mode,'2c') 
          imagesc(handles.out1)
          axis on, colormap(handles.colorcm)
      end 
    end
     
   % set(gcf,'Position',[222    64   944   642])  
      
 
 
    %zlim([0 max(max(handles.stack1(:,:,2)))])
elseif handles.sc ~= 0 %&& handles.ch==0
     imshow(handles.image1)
     caxis([0.5 1])
     %axis on, axis image, colormap(handles.colorcm)
     axis on, colormap(handles.colorcm)
% elseif handles.ch==1 && strcmp(handles.colorcm,'rg')
%  imoverlayNew(double(handles.image1),double(handles.image2),'rg')  
% elseif handles.ch==1 && strcmp(handles.colorcm,'gb')
%  imoverlayNew(double(handles.image1),double(handles.image2),'gb')
%  elseif handles.ch==1 && strcmp(handles.colorcm,'rb')
%  imoverlayNew(double(handles.image1),double(handles.image2),'rb')
end



% --- Executes on selection change in color_choice.
function color_choice_Callback(hObject, eventdata, handles)
% hObject    handle to color_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns color_choice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from color_choice
% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val};
case 'Gray'
   handles.colorcm = 'gray';
case 'Hot'
   handles.colorcm = 'hot';
case 'HSV'
   handles.colorcm = 'hsv';
case 'Jet'
   handles.colorcm = 'jet';
case 'Pink'
   handles.colorcm = 'pink';
case 'Bone'
   handles.colorcm = 'bone';
case '2Ch=rg'
   handles.colorcm = 'rg';
case '2Ch=rb'
   handles.colorcm = 'rb';
case '2Ch=gb'
   handles.colorcm = 'gb';
end
showim(handles)
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function color_choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to color_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in playmov.
function playmov_Callback(hObject, eventdata, handles)
% hObject    handle to playmov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for it = 1:handles.totim
    handles.image1 = handles.stack1(:,:,it);
    if handles.ch==1 && (strcmp(handles.mode,'c') || strcmp(handles.mode,'s'))
     handles.image2 = handles.stack2(:,:,it);
    elseif handles.ch==1 && strcmp(handles.mode,'2c')
     handles.out1=handles.out(:,:,:,it);
    end
    pause(1/handles.totim)
    showim(handles)
    
    %set(gcf,'Position',[222    64   944   642])  
    set(handles.imnb, 'String', it);
    set(handles.imslider,'Value',it/handles.totim)
end
handles.image1 = handles.stack1(:,:,1);
if handles.ch==1;
handles.image2 = handles.stack2(:,:,2);
end
set(handles.imnb, 'String', 1);
set(handles.imslider,'Value',0)
set(handles.imslider,'SliderStep',[handles.playspeed*0.01/handles.totim handles.playspeed*0.1/handles.totim])

showim(handles)


% normalisation function
function matricenorm = normalisation(matrice)
maxi = max(max(matrice));
mini = min(min(matrice));
matricenorm = (matrice-mini)/(maxi-mini);