function varargout = IPM_TnT_Ana(varargin)
% IPM_TNT_ANA MATLAB code for IPM_TnT_Ana.fig
%   IPM_TNT_ANA is a gui script for the analysis of tif stacks derived
%   from kinetic disk images recorded by an andor spinning disk setup
%   utilizing only one camera, which gives a split screen information.
%   It is specifically written as a first step in the analysis of vesicle
%   fusion events using a beam splitter before the camera to split
%   the fluorescence signal into a membrane and a vesicle content channel.
%
%   Furthermore ROIs, which were set utilizing ImageJ can be loaded and
%   coordinates of the ROIs are transferred into both image parts, 
%   allowing to read out the vesicle channel and the occording position
%   in the membrane channel.
%
%   This GUI is not of general interest and is specifically designed
%   to support research within the CRC 803. It was written using MATLab
%   2017b
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IPM_TnT_Ana

% Last Modified by GUIDE v2.5 05-Feb-2019 14:37:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IPM_TnT_Ana_OpeningFcn, ...
                   'gui_OutputFcn',  @IPM_TnT_Ana_OutputFcn, ...
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


% --- Executes just before IPM_TnT_Ana is made visible.
function IPM_TnT_Ana_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for IPM_TnT_Ana
handles.output = hObject;

% Setting up some things withing the GUI
h = gcf;
set(h,'MenuBar','None');
ax = gca;
fig = ancestor(ax, 'figure');
set(ax,'YTick',[]);
set(ax,'XTick',[]);
warning('off');
guidata(hObject, handles);

% UIWAIT makes IPM_TnT_Ana wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IPM_TnT_Ana_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in loadImage.
function loadImage_Callback(hObject, eventdata, handles)
% load a whole time series, which typically consists of several tiff-stacks
% with data sizes arounf 2GB for each stack. To reduce amount of memory
% used information about the whole stack is loaded, specifically which
% frame/time stamp is located in which stack-file, to laod specific images
% only. Reduces amoutn of memory used sramatically, but takes time.
%
% try - catch used only for batch processing, when errors occur these are
% printed for later data evaluation
try 
    % files are loaded and informations are stored to the handles structure
    % of the gui, so the information can be accessed from every callback
    [fn,pn]=uigetfile('*.tif*');
    cd(pn);
    handles.fn=fn;
    handles.pn=pn;
    handles.cposl=[];
    handles.cposr=[];  
    handles.imgStackSlides=[];
    handles.slideFile=[];
    handles.slideNo=[];
    num_images=0;
    fn_base=fn(1:end-5);
    handles.fileNames=[];
    % Inforamtion for every file belonging to this specific timeseries is
    % gathered
    for i=1:str2num(handles.noFiles.String)
        fn=[fn_base,num2str(i-1),'.tif'];
        InfoImage = imfinfo(fn);
        num_images=num_images+numel(InfoImage); % Calcualtion of the total 
                                                % Number of images of this
                                                % timeseries
        handles.imgStackSlides=[handles.imgStackSlides;...
            numel(InfoImage)+sum(handles.imgStackSlides)]; 
        % Storing Information which frame is in which file
        for k=1:numel(InfoImage)
            handles.slideFile=[handles.slideFile;fn];
            handles.slideNo=[handles.slideNo;k];
        end
        handles.fileNames=[handles.fileNames;fn];
    end
    % summing up all information and setting GUI elements like sliders
    % using the total number of frames and such things
    handles.totalImages=num_images;
    handles.slides=[1:num_images];    
    handles.im = imread(fn, 1);
    set(handles.frame,'min', 1);
    set(handles.frame,'max', num_images);
    set(handles.frame,'Value', 1);
    set(handles.frame, 'SliderStep', [1/num_images , 10/num_images]);
    set(handles.systemMessages,'String','Image successfully loaded');
    
    
% Assignments used to verify processes of the gui, uncomment to
% access teh variables in matlab workspace
%     assignin('base','im',handles.im)
%     assignin('base','imslides',handles.imgStackSlides)
%     assignin('base','slides',handles.slides)
%     assignin('base','slideFile',handles.slideFile)
%     assignin('base','fileNames',handles.fileNames)

    % Plotting the first image to fill the gui, update handles
    axes(handles.img);
    ax = gca;
    fig = ancestor(ax, 'figure');
    %imagesc(handles.img,handles.im(:,:,1));
    frame_Callback(hObject, eventdata, handles); % used to plot the image 
                                                 % according to slider
                                                 % position
    axis image;
    set(ax,'YTick',[]);
    set(ax,'XTick',[]);   
catch e
   set(handles.systemMessages,'String','No image is loaded');
   fprintf(1,'The identifier was:\n%s',e.identifier);
   fprintf(1,'There was an error! The message was:\n%s',e.message);
end
guidata(hObject, handles);

% --- Executes on slider movement.
function frame_Callback(hObject, eventdata, handles)
% used for plotting, when the slider is moved, also called to plot actual
% image by other functions
global tpos %stores position, is decapricated as slider info is now used
 try
     % Setting up the plot
    axes(handles.img)
    ax = gca;
    axis('manual');
    frameIdx=floor(get(handles.frame,'Value')); % Frame to plot by slider
                                                % position.
    slIdx=find(handles.slides==frameIdx);
    handles.slideFile(frameIdx,:)
    % Reading stack and specific image using data stored before
    handles.im=imread(handles.slideFile(frameIdx,:),...
        handles.slideNo(frameIdx));
    % Checking if user wants to see the total sum of the image stack, if
    % stack is summed up and shown.
    if handles.stackSum.Value==1
       Slices=str2num(handles.NSlice.String);
       frameIdx=floor(get(handles.frame,'Value'))
       StartSlice=frameIdx-Slices;
       if StartSlice<1
           StartSlice=1;
       end
       EndSlice=frameIdx+Slices;
       % Just to be sure when slider gives a slide after the max slide
       % number
       if EndSlice>=max(handles.slides)
           EndSlice=max(handles.slides);
       end
       % Checking if slides to be summed are spanned across several files       
       if EndSlice-StartSlice>3000 % derived from the standard values of 
                                   % the andor software
           for o=1:size(handles.fileNames,1)
               % Reading the images for sum calculation
                FileTif=handles.fileNames(o,:);
                InfoImage=imfinfo(FileTif);
                mImage=InfoImage(1).Width;
                nImage=InfoImage(1).Height;
                NumberImages=length(InfoImage);
                handles.im=zeros(nImage,mImage,NumberImages,'uint16');
                TifLink = Tiff(FileTif, 'r');
                for i=1:NumberImages
                   TifLink.setDirectory(i);
                   dmy(:,:,i)=TifLink.read();
                end
                dmy2(:,:,o)=sum(dmy,3);
                %assignin('base','dmy2',dmy2) % to check just in case
                handles.im=sum(dmy2,3);
                TifLink.close();
           end
           clear('dmy','dmy2');
       else
           im=uint16([zeros(size(handles.im))]);
           t=1;
           for p=StartSlice:EndSlice           
                im(:,:,t)=imread(handles.slideFile(p,:),handles.slideNo(p));
                t=t+1;
           end
           %assignin('base','im',im) % to check just in case
           handles.im=sum(im,3);   
           %assignin('base','him',handles.im) % to check just in case
       end
    end   
    
    % Plotting the resulst, if sum chosen, the summed stacked; if not image
    % according to slider position. Also user colormap is chosen according
    % to GUI settings by user.
    imagesc(handles.img,handles.im);
    c=cellstr(get(handles.setColorMap,'String'));
    cidx=get(handles.setColorMap,'Value');
    colormap(c{cidx});   
    axis image;
    set(ax,'YTick',[]);
    set(ax,'XTick',[]);
    hold on
    axis('manual');
    frameStamps={};
    drawnow
    
    % When user wants to see a time stamp/ frame number
    if handles.showTimestamp.Value==1
        axes(handles.graph2);
        ax=gca;
        fig = ancestor(ax, 'figure');   
        %plot(handles.Traces_r(:,tpos),'r-');        
        frameIdx=floor(get(handles.frame,'Value'))
        hold on
        vline(frameIdx,'k')
        hold off   
        drawnow        
        axes(handles.graph);        
        ax=gca;        
        fig = ancestor(ax, 'figure');   
        %plot(handles.Traces_l(:,tpos),'b-');
        frameIdx=floor(get(handles.frame,'Value'));     
        hold on        
        vline(frameIdx,'k')
        hold off
        drawnow       
    end    
    % as always update the handles, so every function has fresh
    % information!
    guidata(hObject, handles);
 catch
    set(handles.systemMessages,'String','Something went wront while plotting the frame')
 end
hold off

% --- Executes during object creation, after setting all properties.
function frame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function micrometerPerPixel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function micrometerPerPixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in setColorMap.
function setColorMap_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
frame_Callback(hObject, eventdata, handles) %just plotting when colormap is
                                            %changed
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function setColorMap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Deprecated as users found to be unnecessary
% % --- Executes on button press in saveMovie.
% function saveMovie_Callback(hObject, eventdata, handles)
% % hObject    handle to saveMovie (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% try
%     set(handles.systemMessages,'String','Saving movie in progress')
%     F=[];
%     fn=get(handles.filenameSave,'String');
%     if get(handles.timeStamp,'Value')==1
%         c=clock;
%         fn_mov=strcat(fn,'_',num2str(c(1)),'_',num2str(c(2)),'_',...
%             num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)),'_',...
%             num2str(c(6)),'.avi');
%     else
%         fn_mov=strcat(fn,'avi');
%     end
%     v = VideoWriter(fn_mov,'Uncompressed AVI');
%     v.FrameRate=10;
%     %v.Quality=100;    
%     open(v)
%     for i=1:length(handles.im)
%         F = [F;getframe(handles.img)];
%         set(handles.frame,'Value',i);
%         guidata(hObject, handles);
%         frame_Callback(hObject, eventdata, handles);        
%         guidata(hObject, handles);        
%     end 
%     writeVideo(v,F)
%     close(v)
%     %movie2avi(F,'image.avi','Compression','Cinepak')
%     guidata(hObject, handles);
%     set(handles.systemMessages,'String','Movie saved')
% catch
%    set(handles.systemMessages,'String','Was not able to save movie')
% end

% Deprecated as found to be unnecessary by users
% --- Executes on selection change in loadMethod.
function loadMethod_Callback(hObject, eventdata, handles)
g=cellstr(get(handles.loadMethod,'String'));
gidx=get(handles.loadMethod,'Value');
if strcmp(g{gidx},'large')
    set(handles.systemMessages,'String','Loading one slice at a time, might be very slow');
end

% --- Executes during object creation, after setting all properties.
function loadMethod_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stackSum.
function stackSum_Callback(hObject, eventdata, handles)
guidata(hObject,handles)
frame_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

function NSlice_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function NSlice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alignCoords.
% Copy Coordinates to second channel, Instructions are given within the GUI
function alignCoords_Callback(hObject, eventdata, handles)
set(handles.systemMessages,'String','Use left mouse button to altering points in both image sides, First left, then right!')
try       
    % first clean up matrices with centerpositions
    handles.cposl=rmmissing(handles.alignment.Data(:,1:2));
    handles.cposr=rmmissing(handles.alignment.Data(:,3:4));
    % when center positions are already there, make cells and empty them
    % could be done much easier! This justr cleans up when user does it a
    % second time
    if iscell(handles.cposl)
        handles.cposl=cell2mat(handles.cposl);
        handles.cposr=cell2mat(handles.cposr);
        handles.cposl=[];
        handles.cposr=[];
    end    
    % setting up coordinate copy from one window part to the other using
    % features visible in both channels. Users clicks on these features
    % and the coordinate shift is calculated. Gets more exact when more
    % features are used
    ax = gca;
    fig = ancestor(ax, 'figure');
    axis('manual');
    hold on
    for m=1:2
        c = ginput(1);
        sel = get(fig, 'SelectionType');
        scatter(c(1),c(2),'wx');
        if m==1
            handles.cposl=[handles.cposl;c];        
        else
            handles.cposr=[handles.cposr;c];        
        end
    end
    % Store results of ROI alignment for both channels
    handles.alignment.Data=[handles.cposl,handles.cposr];    
    hold off
    hold on
    % show what was done
    scatter(handles.cposl(:,1),handles.cposl(:,2),'wx');
    scatter(handles.cposr(:,1),handles.cposr(:,2),'yx');
    hold off
    % calculating the coordinate shofts from every position chosen and 
    % claculating the mean 
    diffPos=handles.cposr-handles.cposl;
    handles.Xdiffmean=-mean(diffPos(:,1));
    handles.Ydiffmean=-mean(diffPos(:,2));
catch err          
    set(handles.systemMessages,'String','Something unexpected happened, try setting centers again. Follow the instructions!');
    warndlg(err.message);
    disp(err)
end
guidata(hObject, handles);


% --- Executes on button press in setRoi.
function setRoi_Callback(hObject, eventdata, handles)
% set ROIs on one channel, which are copied to the second channel
% only necessary if user after IMAGEJ wants to add some ROIs.
    handles.rpos=rmmissing(handles.roiTable.Data(:,1:2));    
    if iscell(handles.rpos)
        handles.rpos=cell2mat(handles.rpos);
        handles.rpos=[];
    end    
    ax = gca;
    fig = ancestor(ax, 'figure');
    axis('manual');
    hold on
    r = ginput(1);
    sel = get(fig, 'SelectionType');
    scatter(r(1),r(2),'rx');
    scatter(r(1)+handles.Xdiffmean,r(2)+handles.Ydiffmean,'mx');
    % Roi position and radius are stored
    handles.rpos=[handles.rpos;r];
    handles.roiTable.Data=[handles.rpos,ones(size(handles.rpos,1)).*3];       
    hold off
    hold on
    scatter(handles.rpos(:,1),handles.rpos(:,2),'r.');
    hold off


% --- Executes on button press in showRoi.
function showRoi_Callback(hObject, eventdata, handles)
% show all ROIs!
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)
frame_Callback(hObject, eventdata, handles);
for i=1:size(handles.roiTable.Data,1)
    ang=0:0.01:2*pi; 
    xp=handles.roiTable.Data(i,3).*cos(ang);
    yp=handles.roiTable.Data(i,3).*sin(ang);
    hold on
    plot(handles.roiTable.Data(i,1)+xp,handles.roiTable.Data(i,2)+yp,'w-')
    plot(handles.roiTable.Data(i,4)+xp,handles.roiTable.Data(i,5)+yp,'w-');
    %plot(handles.Xdiffmean+handles.roiTable.Data(i,1)+xp,...
    %    handles.Ydiffmean+handles.roiTable.Data(i,2)+yp,'r-');
    hold off
end


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% Saving results to results.mat for further analysis using Fusionanalyzer
results=[];
results.Trace_BG_l=handles.Traces_BG_l;
results.Trace_l=handles.Traces_l;
results.Trace_BG_r=handles.Traces_BG_r;
results.Trace_r=handles.Traces_r;
results.ROIs=handles.roiTable.Data(:,1:5);
results.Alignment=handles.alignment.Data;
results.XdiffROI=handles.Xdiffmean;
results.YdiffROI=handles.Ydiffmean;
results.fileNames=handles.fileNames;
save('results.mat','results')


% --- Executes on button press in playMovie.
function playMovie_Callback(hObject, eventdata, handles)
% Sometime you just want to see whats happening in the time series
global stop
stop = 0;
for i=get(handles.frame,'Value'):max(handles.slides)   
    set(handles.frame,'Value',i);
    guidata(hObject, handles);
    frame_Callback(hObject, eventdata, handles);        
    guidata(hObject, handles); 
    drawnow
    if stop==1
        break
    end
end


% --- Executes on mouse press over axes background.
function img_ButtonDownFcn(hObject, eventdata, handles)
% Checks for click in image while movie is running, to stop the movie
global stop
stop=1


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% just to stop the movie while running
global stop
stop=1


function noFiles_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function noFiles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadImageJROI.
function loadImageJROI_Callback(hObject, eventdata, handles)
% Loading ROI files from IMAGEJ and storing these in the handles.
% typically these are zip-files, these are unzipped in the current folder
% Therefore this GUI- files should be located in the same folder where the
% data is located. As the PhD typically make changes to the GUI behaviour
% to ease data analysis this is normally done. For a general Program one
% should use defined paths
[fn,pn]=uigetfile('*.zip');
mkdir('.\ImageJROI\');
unzip(strcat(pn,fn),'.\ImageJROI\');
ROIlist=dir('.\ImageJROI\');
assignin('base','roilist',ROIlist)
handles.roiTable.Data=[];
guidata(hObject, handles);
for i=3:size(ROIlist,1)
    roi=ReadImageJROI(strcat('.\ImageJROI\',ROIlist(i).name));
    y=(roi.vnRectBounds(1)+roi.vnRectBounds(3))./2;
    x=(roi.vnRectBounds(2)+roi.vnRectBounds(4))./2;
    r=abs((roi.vnRectBounds(1)-roi.vnRectBounds(3))./2);
    handles.roiTable.Data=[handles.roiTable.Data;[x y r]];
    guidata(hObject, handles);
end
roiTrans_Callback(hObject, eventdata, handles)
guidata(hObject, handles);




% --- Executes on button press in contrast.
function contrast_Callback(hObject, eventdata, handles)
% To enhance contrast in the images, if nothing can be seen, typically 
% not used
guidata(hObject, handles);
imcontrast(handles.img)
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of contrast


% --- Executes on button press in readRoi.
% Reading imageJ stored ROIs Intensities
function readRoi_Callback(hObject, eventdata, handles)

intensity_l=[];
intensity_r=[];
results={};

ROI_coords=[];
ROI_coords=handles.roiTable.Data(:,1:5);
%ROI_coords(:,4)=ROI_coords(:,1)+handles.Xdiffmean;
%ROI_coords(:,5)=ROI_coords(:,2)+handles.Ydiffmean;    

% Creating masks from the ROI information and reading the intensities of
% both channels, calculating the mean within the ROI
tic
for i=1:size(handles.fileNames,1) 
    FileTif=handles.fileNames(i,:);
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    im=zeros(nImage,mImage,NumberImages,'uint16');
    TifLink = Tiff(FileTif, 'r');
    for k=1:NumberImages
       TifLink.setDirectory(k);
       im(:,:,k)=TifLink.read();
    end
    TifLink.close();
        
    for t=1:size(handles.roiTable.Data(:,1:3),1)                
        ROI_masks_l=createCirclesMask(handles.im,ROI_coords(t,4:5),ROI_coords(t,3));
        ROI_masks_r=createCirclesMask(handles.im,ROI_coords(t,1:2),ROI_coords(t,3));
        for p=1:NumberImages
            % Calculating the mean value of the ROI
            intensity_l=[intensity_l;mean(mean(double(im(:,:,p)).*double(ROI_masks_l)))];        
            intensity_r=[intensity_r;mean(mean(double(im(:,:,p)).*double(ROI_masks_r)))];        
        end
        results{t,i,1}=intensity_l;
        results{t,i,2}=intensity_r;    
        intensity_r=[];
        intensity_l=[];
    end    
end
% Using the mean of every ROI in every frame to create an intensitiy trace
handles.Traces_r=[];
handles.Traces_l=[];
for e=1:size(results,1)
    trace_r=[];
    trace_l=[];
    for f=1:size(results,2)
        trace_r=[trace_r;results{e,f,2}];
        trace_l=[trace_l;results{e,f,1}];
    end    
    handles.Traces_r=[handles.Traces_r,trace_r];
    handles.Traces_l=[handles.Traces_l,trace_l];
end
toc
guidata(hObject, handles);
axes(handles.graph)
ax = gca;
axis('manual');
% All Traces shown are background corrected, but only when shown
if handles.bgRight.Value==1
    plot(handles.Traces_r(:,1)-handles.Traces_BG_r,'r-');
else
    plot(handles.Traces_r(:,1),'r-');
end
drawnow
guidata(hObject, handles);
axes(handles.graph2)
ax = gca;
axis('manual');
if handles.bgLeft.Value==1
    plot(handles.Traces_l(:,1)-handles.Traces_BG_l,'bo');
else
    plot(handles.Traces_l(:,1),'bo');
end

drawnow
assignin('base','results',results)
assignin('base','Traces_r',handles.Traces_r)
guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in roiTable.
function roiTable_CellSelectionCallback(hObject, eventdata, handles)
% Table of ROIs read from ImageJ shown in the GUI and possible
% manipulation of the ROIs by the user, like deleting one from the table in
% the GUI. ROIs are replotted on the image when something in the table is
% changed
global tpos
frame_Callback(hObject, eventdata, handles);
handles.roiTable.Data(eventdata.Indices(1),1:3)
tpos=eventdata.Indices(1);
guidata(hObject, handles);
axes(handles.graph);
ax = gca;
hold off
%axis('manual');
if handles.bgRight.Value==0
    plot(handles.Traces_r(:,eventdata.Indices(1)),'b.');
else
    plot(handles.Traces_r(:,eventdata.Indices(1))-handles.Traces_BG_r,'b.');
end
drawnow
guidata(hObject, handles);
axes(handles.graph2);
ax = gca;
hold off
if handles.bgLeft.Value==0    
    plot(handles.Traces_l(:,eventdata.Indices(1)),'r.');
else
    plot(handles.Traces_l(:,eventdata.Indices(1))-handles.Traces_BG_l,'r.');
end
drawnow
axes(handles.img);
ax = gca;
hold on
ang=0:0.01:2*pi; 
xp=handles.roiTable.Data(eventdata.Indices(1),3).*cos(ang);
yp=handles.roiTable.Data(eventdata.Indices(1),3).*sin(ang);
hold on
plot(handles.roiTable.Data(eventdata.Indices(1),1)+xp,handles.roiTable.Data(eventdata.Indices(1),2)+yp,'m-')
plot(handles.roiTable.Data(eventdata.Indices(1),4)+xp,handles.roiTable.Data(eventdata.Indices(1),5)+yp,'m-');
%plot(handles.Xdiffmean+handles.roiTable.Data(i,1)+xp,...
%    handles.Ydiffmean+handles.roiTable.Data(i,2)+yp,'r-');
    hold off
    if handles.showTimestamp.Value==1
        axes(handles.graph);
        ax=gca;
        fig = ancestor(ax, 'figure');   
        cla(ax)
        drawnow
        if handles.bgRight.Value==0
            plot(handles.Traces_r(:,eventdata.Indices(1)),'b.');
        else
            plot(handles.Traces_r(:,eventdata.Indices(1))-handles.Traces_BG_r,'b.');
        end       
        frameIdx=floor(get(handles.frame,'Value'))
        hold on
        vline(frameIdx,'k')
        hold off   
        drawnow     
        
        axes(handles.graph2);        
        ax=gca;        
        fig = ancestor(ax, 'figure');   
        cla(ax)
        drawnow
        if handles.bgLeft.Value==0    
            plot(handles.Traces_l(:,eventdata.Indices(1)),'r.');
        else
            plot(handles.Traces_l(:,eventdata.Indices(1))-handles.Traces_BG_l,'r.');
        end
            frameIdx=floor(get(handles.frame,'Value'));     
        hold on        
        vline(frameIdx,'k')
        hold off
        drawnow       
    end
guidata(hObject, handles);


% --- Executes on button press in showTimestamp.
function showTimestamp_Callback(hObject, eventdata, handles)
    guidata(hObject,handles)
    frame_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of showTimestamp


% --- Executes on selection change in roiTranslation.
function roiTranslation_Callback(hObject, eventdata, handles)
% Just Press to copy the ROIs from one half image to the other half, to
% enable position readout for both channels
    diffPos=handles.cposr-handles.cposl;
    
    c=cellstr(get(handles.roiTranslation,'String'));
    cidx=get(handles.roiTranslation,'Value');
    switch c{cidx}
        case 'left'
            handles.Xdiffmean=-mean(diffPos(:,1));
            handles.Ydiffmean=-mean(diffPos(:,2));
        case 'right'
            handles.Xdiffmean=mean(diffPos(:,1));
            handles.Ydiffmean=mean(diffPos(:,2));
    end       

    
% Hints: contents = cellstr(get(hObject,'String')) returns roiTranslation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiTranslation


% --- Executes during object creation, after setting all properties.
function roiTranslation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiTranslation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roiTrans.
function roiTrans_Callback(hObject, eventdata, handles)
% hObject    handle to roiTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roiTable.Data(:,4)=handles.roiTable.Data(:,1)-handles.Xdiffmean;
handles.roiTable.Data(:,5)=handles.roiTable.Data(:,2)-handles.Ydiffmean;
guidata(hObject, handles);
showRoi_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes on button press in setBGRoi.
function setBGRoi_Callback(hObject, eventdata, handles)
% Setting and reading out ROI for Background correction, translation into
% both channels is done like for the ReadOut ROIs
    handles.BGpos_l=[];
    handles.BGpos_r=[];
    ax = gca;
    fig = ancestor(ax, 'figure');
    axis('manual');
    hold on
    r = ginput(1);
    sel = get(fig, 'SelectionType');
    scatter(r(1),r(2),'yo');
    scatter(r(1)-handles.Xdiffmean,r(2)-handles.Ydiffmean,'yo');
    handles.BGpos_r=r;
    handles.BGpos_l=[r(1)-handles.Xdiffmean,r(2)-handles.Ydiffmean];    
    hold off
    
    intensity_l=[];
    intensity_r=[];
    handles.Traces_BG_r=[];
    handles.Traces_BG_l=[];
    
for i=1:size(handles.fileNames,1) 
    FileTif=handles.fileNames(i,:);
    InfoImage=imfinfo(FileTif);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    im=zeros(nImage,mImage,NumberImages,'uint16');
    TifLink = Tiff(FileTif, 'r');
    for k=1:NumberImages
       TifLink.setDirectory(k);
       im(:,:,k)=TifLink.read();
    end
    TifLink.close();       
            
        ROI_masks_bg_l=createCirclesMask(handles.im,handles.BGpos_l,2);
        ROI_masks_bg_r=createCirclesMask(handles.im,handles.BGpos_r,2);
        for p=1:NumberImages
            intensity_l=[intensity_l;mean(mean(double(im(:,:,p)).*double(ROI_masks_bg_l)))];        
            intensity_r=[intensity_r;mean(mean(double(im(:,:,p)).*double(ROI_masks_bg_r)))];        
        end
        results{i,1}=intensity_l;
        results{i,2}=intensity_r;    
        intensity_r=[];
        intensity_l=[];
  
end

handles.BG_Traces_r=[];
handles.BG_Traces_l=[];

trace_r=[];
trace_l=[];
for f=1:size(results,1)
    trace_r=[trace_r;results{f,2}];
    trace_l=[trace_l;results{f,1}];
end
handles.Traces_BG_r=[handles.Traces_BG_r,trace_r];
handles.Traces_BG_l=[handles.Traces_BG_l,trace_l];
assignin('base','bg_l',handles.Traces_BG_l)
assignin('base','bg_r',handles.Traces_BG_r)
guidata(hObject, handles);


% --- Executes on button press in bgLeft.
function bgLeft_Callback(hObject, eventdata, handles)
% hObject    handle to bgLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgLeft


% --- Executes on button press in bgRight.
function bgRight_Callback(hObject, eventdata, handles)
% hObject    handle to bgRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgRight
