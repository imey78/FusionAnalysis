function varargout = FusionAnalyser(varargin)
% FUSIONANALYSER MATLAB code for FusionAnalyser.fig
%   Fusionanalyzer is a gui script for the analysis of time traces derived
%   by IPM_TNT_Ana. This Gui was specifically designed to evaluate on data 
%   and experiments derived in the Lab auf Prof. Steinem. It is meant to
%   deal with a large data set, stored in specific locations and is not of
%   general interest. This GUI derives data from a specific folder
%   structure and allow the scientist to create a "judgment" string, which
%   categorizes single vesicle fusion events on pore spanning membranes.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FusionAnalyser

% Last Modified by GUIDE v2.5 18-Mar-2019 14:04:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FusionAnalyser_OpeningFcn, ...
                   'gui_OutputFcn',  @FusionAnalyser_OutputFcn, ...
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


% --- Executes just before FusionAnalyser is made visible.
function FusionAnalyser_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes FusionAnalyser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FusionAnalyser_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% loading a specific time series, which might be splitted into several
% files and the according results.mat of IPM_TNT_ANA
% the time series is plotted as well as the intensity of the rois for these
% time series.
f_img=uigetfile('*.tif');
f_res=uigetfile('*.mat');
fn_base=f_img(1:end-5);
handles.fileNames=[];
info=imfinfo([fn_base,'1.tif']); % using secong stack, but could be every 
                                 % one. All time series measured by Peter
                                 % used to be exactly 6 stacks 
dmy=zeros(info.Width,info.Height);
for i=1:str2num(handles.lfNo.String)
    fn=[fn_base,num2str(i-1),'.tif'];
    dmy=dmy+imread(fn); %unnecessary complicated but works
end
axes(handles.axes1);
ax = gca;
imagesc(dmy) % plot the image
load(f_res); % loading according results.mat
handles.ROIs=results.ROIs; % getting the set rois into a GUI table
handles.roiTable.Data=handles.ROIs;
handles.judge_string=cell(size(handles.ROIs,1),1) % create the string which 
                                                  % will descirbe which
                                                  % type of event the user
                                                  % looking at
handles.judge_string
% loading all traces of results.mat into GUI handles
handles.traces_l=results.Trace_l;
handles.traces_r=results.Trace_r;
handles.bg_l=results.Trace_BG_l;
handles.bg_r=results.Trace_BG_r;
guidata(hObject, handles);

function lfNo_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function lfNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showROI.
function showROI_Callback(hObject, eventdata, handles)
% Just when the user wants to see which ROIs were set
for i=1:size(handles.roiTable.Data,1)
    ang=0:0.01:2*pi; 
    xp=handles.roiTable.Data(i,3).*cos(ang);
    yp=handles.roiTable.Data(i,3).*sin(ang);
    hold on
    plot(handles.roiTable.Data(i,1)+xp,handles.roiTable.Data(i,2)+yp,'m-')
    plot(handles.roiTable.Data(i,4)+xp,handles.roiTable.Data(i,5)+yp,'m-');
    %plot(handles.Xdiffmean+handles.roiTable.Data(i,1)+xp,...
    %    handles.Ydiffmean+handles.roiTable.Data(i,2)+yp,'r-');
    hold off
end
guidata(hObject, handles);

% --- Executes when selected cell(s) is changed in roiTable.
function roiTable_CellSelectionCallback(hObject, eventdata, handles)
% When ROIs are selected within the GUI Table, the ROI is plotted and the
% according traces are shown
global tpos
axes(handles.axes1);
ax = gca;
for i=1:size(handles.roiTable.Data,1)
    ang=0:0.01:2*pi; 
    xp=handles.roiTable.Data(i,3).*cos(ang);
    yp=handles.roiTable.Data(i,3).*sin(ang);
    hold on
    plot(handles.roiTable.Data(i,1)+xp,handles.roiTable.Data(i,2)+yp,'m-')
    plot(handles.roiTable.Data(i,4)+xp,handles.roiTable.Data(i,5)+yp,'m-');
    %plot(handles.Xdiffmean+handles.roiTable.Data(i,1)+xp,...
    %    handles.Ydiffmean+handles.roiTable.Data(i,2)+yp,'r-');
    hold off
end
tpos=eventdata.Indices(1);
handles.roiTable.Data(eventdata.Indices(1),1:3)
hold on
ang=0:0.01:2*pi; 
xp=handles.roiTable.Data(eventdata.Indices(1),3).*cos(ang);
yp=handles.roiTable.Data(eventdata.Indices(1),3).*sin(ang);
hold on
plot(handles.roiTable.Data(eventdata.Indices(1),1)+xp,handles.roiTable.Data(eventdata.Indices(1),2)+yp,'k-')
plot(handles.roiTable.Data(eventdata.Indices(1),4)+xp,handles.roiTable.Data(eventdata.Indices(1),5)+yp,'k-');
%plot(handles.Xdiffmean+handles.roiTable.Data(i,1)+xp,...
%    handles.Ydiffmean+handles.roiTable.Data(i,2)+yp,'r-');
hold off

%Traces shown are corrected by the background fluorescence
axes(handles.axes2);
ax = gca;
if handles.BG_Corr_l.Value==1
    plot(handles.traces_l(:,eventdata.Indices(1))-handles.bg_l)
    trace_l=(handles.traces_l(:,eventdata.Indices(1))-handles.bg_l);
else
    plot(handles.traces_l(:,eventdata.Indices(1)))
    trace_l=(handles.traces_l(:,eventdata.Indices(1)));
end

axes(handles.axes3);
ax = gca;
if handles.BG_Corr_r.Value==1
    plot(handles.traces_r(:,eventdata.Indices(1))-handles.bg_r)
else
    plot(handles.traces_r(:,eventdata.Indices(1)))
end

axes(handles.axes2);
ax = gca;
% To cope for crosstalk between channels a crosstalk correction is done,
% further information in Publication SI
cf=str2num(handles.ct_factor.String)
handles.crosstalk.Value
if handles.crosstalk.Value==1
    handles.BG_Corr_l.Value
    if ((handles.BG_Corr_r.Value==1) && (handles.BG_Corr_l.Value==1))
        disp('test')
        plot((handles.traces_l(:,eventdata.Indices(1))-handles.bg_l)-((handles.traces_r(:,eventdata.Indices(1))-handles.bg_r).*cf))
        trace_l=((handles.traces_l(:,eventdata.Indices(1))-handles.bg_l)-((handles.traces_r(:,eventdata.Indices(1))-handles.bg_r).*cf));
    elseif handles.BG_Corr_r.Value==0 && handles.BG_Corr_l.Value==0
        plot((handles.traces_l(:,eventdata.Indices(1)))-((handles.traces_r(:,eventdata.Indices(1))).*cf))
        trace_l=((handles.traces_l(:,eventdata.Indices(1)))-((handles.traces_r(:,eventdata.Indices(1))).*cf));
    elseif handles.BG_Corr_r.Value==1 && handles.BG_Corr_l.Value==0
        plot((handles.traces_l(:,eventdata.Indices(1)))-((handles.traces_r(:,eventdata.Indices(1))-handles.bg_r).*cf))
        trace_l=((handles.traces_l(:,eventdata.Indices(1)))-((handles.traces_r(:,eventdata.Indices(1))-handles.bg_r).*cf));
    elseif handles.BG_Corr_r.Value==0 && handles.BG_Corr_l.Value==1    
        plot((handles.traces_l(:,eventdata.Indices(1))-handles.bg_l)-((handles.traces_r(:,eventdata.Indices(1))).*cf))
        trace_l=((handles.traces_l(:,eventdata.Indices(1))-handles.bg_l)-((handles.traces_r(:,eventdata.Indices(1))).*cf));
    end
end

% Several options to filter teh shown Signals to ease identification of
% event type
switch handles.filterChoice.String{handles.filterChoice.Value}
    case 'none'
%         axes(handles.axes2);
%         ax = gca;
%         if handles.BG_Corr_l.Value==1
%             plot(handles.traces_l(:,tpos)-handles.bg_l)
%         else
%             plot(handles.traces_l(:,tpos))
%         end
%         
%         
% 
%         axes(handles.axes3);
%         ax = gca;
%         if handles.BG_Corr_r.Value==1
%             plot(handles.traces_r(:,tpos)-handles.bg_r)
%         else
%             plot(handles.traces_r(:,tpos))
%         end
%         guidata(hObject, handles);
    case 'mean'        
        axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            dmy=movmean(trace_l,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmean(trace_l,str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            dmy=movmean(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmean(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off
    case 'median'
        axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            dmy=movmedian(trace_l,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmedian(trace_l,str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            dmy=movmedian(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmedian(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off
    case 'sq'
    case 'envelope'
                axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            [dmy_u,dmy_l]=envelope(handles.traces_l(:,tpos)-handles.bg_l,str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');
        else
            [dmy_u,dmy_l]=envelope(handles.traces_l(:,tpos),str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            [dmy_u,dmy_l]=envelope(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');            
        else
            [dmy_u,dmy_l]=envelope(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');
        end
        hold off
end
guidata(hObject, handles);


% --- Executes on button press in BG_Corr_l.
function BG_Corr_l_Callback(hObject, eventdata, handles)


% --- Executes on button press in BG_Corr_r.
function BG_Corr_r_Callback(hObject, eventdata, handles)


% --- Executes on selection change in filterChoice.
% Several options to filter teh shown Signals to ease identification of
% event type, all traces shown are corrected by the background ROI
function filterChoice_Callback(hObject, eventdata, handles)
global tpos
axes(handles.axes2);
ax = gca;
if handles.BG_Corr_l.Value==1
    plot(handles.traces_l(:,tpos)-handles.bg_l)
else
    plot(handles.traces_l(:,tpos))
end

axes(handles.axes3);
ax = gca;
if handles.BG_Corr_r.Value==1
    plot(handles.traces_r(:,tpos)-handles.bg_r)
else
    plot(handles.traces_r(:,tpos))
end
guidata(hObject, handles);


switch handles.filterChoice.String{handles.filterChoice.Value}
    case 'none'
        axes(handles.axes2);
        ax = gca;
        if handles.BG_Corr_l.Value==1
            plot(handles.traces_l(:,tpos)-handles.bg_l)
        else
            plot(handles.traces_l(:,tpos))
        end

        axes(handles.axes3);
        ax = gca;
        if handles.BG_Corr_r.Value==1
            plot(handles.traces_r(:,tpos)-handles.bg_r)
        else
            plot(handles.traces_r(:,tpos))
        end
        guidata(hObject, handles);
    case 'mean'        
        axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            dmy=movmean(handles.traces_l(:,tpos)-handles.bg_l,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmean(handles.traces_l(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            dmy=movmean(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmean(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off
    case 'median'
        axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            dmy=movmedian(handles.traces_l(:,tpos)-handles.bg_l,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmedian(handles.traces_l(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            dmy=movmedian(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy,'r-');
        else
            dmy=movmedian(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy,'r-')
        end
        hold off
    case 'sq'
    case 'envelope'
        axes(handles.axes2);
        ax = gca;
        hold on
        if handles.BG_Corr_l.Value==1
            [dmy_u,dmy_l]=envelope(handles.traces_l(:,tpos)-handles.bg_l,str2num(handles.filtSize.String));
            plot(dmy_u,'r-',dmy_l);
        else
            [dmy_u,dmy_l]=envelope(handles.traces_l(:,tpos),str2num(handles.filtSize.String));
            plot(dmy_u,'r-',dmy_l);
        end
        hold off

        axes(handles.axes3);
        ax = gca;
        hold on
        if handles.BG_Corr_r.Value==1
            [dmy_u,dmy_l]=envelope(handles.traces_r(:,tpos)-handles.bg_r,str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');
        else
            [dmy_u,dmy_l]=envelope(handles.traces_r(:,tpos),str2num(handles.filtSize.String));
            plot(dmy_u,'r-');
            plot(dmy_l,'r-');
        end
        hold off
end


% --- Executes during object creation, after setting all properties.
function filterChoice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtSize_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function filtSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in judge.
% First hierarchie of judging the events
function judge_Callback(hObject, eventdata, handles)
global count
count=1;
str= "u: use; " + newline + "r: rubbish; " +newline + "l: later; ";
handles.msg.String=str;
guidata(hObject, handles);
%  switch eventdata.Key
%     case 'u'
%         waitforbuttonpress
%         eventdata.Key
%     case 'r'
%         judge_Callback(handles.judge, [], handles); 
%     case 'r'
%         judge_Callback(handles.judge, [], handles);    
%  end


% --- Executes on key press with focus on figure1 and none of its controls.
% To initiae event categorization press j on figure
function figure1_KeyPressFcn(hObject, eventdata, handles)
% count is the hierarchie of the decision tree underlying the events to be
% categorized in, a decisive string is built up describing the type of
% event shown in the ROI selected. identified by the user looking at the
% possibly filtered traces
global count
global tpos
switch eventdata.Key
    case 'j'
        judge_Callback(handles.judge, [], handles); 
end

if count==1
str= "u: use; " + newline + "r: rubbish; " +newline + "l: later; ";
handles.msg.String=str;
guidata(hObject, handles);
 switch eventdata.Key
    case 'u'
        tpos
        eventdata.Key
        count=2;
        str= "p: fusion pore" + newline + "n: no fusion pore";
        handles.msg.String=str;
        guidata(hObject, handles);
    case 'r'
        handles.judge_string{tpos}='r'            
        guidata(hObject, handles);
        judge_Callback(hObject, eventdata, handles)
    case 'l'
        handles.judge_string{tpos}='l'            
        guidata(hObject, handles);
        judge_Callback(hObject, eventdata, handles)           
 end
elseif count==2
    switch eventdata.Key
    case 'p'
        tpos
        eventdata.Key
        count=3;
        str= "f: full fusion" + newline + "i: incomplete fusion";
        handles.msg.String=str;
        guidata(hObject, handles);
    case 'n'
        tpos
        eventdata.Key
        count=3;
        str= "b: burst" + newline + "e: dead end hemifusion" + newline + "d: pure docking";
        handles.msg.String=str;
        guidata(hObject, handles);
    end
elseif count==3
    switch eventdata.Key
        case 'f'
            count=4;
            str= "d: direct fusion" + newline + "p: sticking pore" + newline + "u: unstable hemifusion diaphragm" + newline + "s: stable hemifusion diaphragm";
            handles.msg.String=str;
            guidata(hObject, handles);
        case 'i'
            count=4;
            str= "c: fusion pore closing" +  newline + "a: vesicle aggregate";
            handles.msg.String=str;
            guidata(hObject, handles);
        case 'b'
            handles.judge_string{tpos}='unb'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'e'
            handles.judge_string{tpos}='une'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'd'
            handles.judge_string{tpos}='und'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
    end
    
elseif count==4
    switch eventdata.Key
        case 'd'
            handles.judge_string{tpos}='upfd'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'p'            
            handles.judge_string{tpos}='upfp'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'u'            
            handles.judge_string{tpos}='upfu'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 's'            
            handles.judge_string{tpos}='upfs'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'c'            
            handles.judge_string{tpos}='upic'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
        case 'a'            
            handles.judge_string{tpos}='upia'            
            guidata(hObject, handles);
            judge_Callback(hObject, eventdata, handles)
    end
end
guidata(hObject, handles);
 assignin('base','judge_string',handles.judge_string)

     
     


% --- Executes on button press in saveJudge.
function saveJudge_Callback(hObject, eventdata, handles)
% saving the event categorization into the judge.mat
jstr=handles.judge_string;
save('judge.mat','jstr')


% --- Executes on button press in crosstalk.
function crosstalk_Callback(hObject, eventdata, handles)


function ct_factor_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ct_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
