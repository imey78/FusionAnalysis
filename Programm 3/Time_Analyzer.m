function varargout = Time_Analyzer(varargin)
% TIME_ANALYZER MATLAB code for Time_Analyzer.fig
%
% This is the code used for the analysis of time differences between the
% two Channel Signals derived by FusionAnalyzer and IPM_ANA_TNT. IT is
% basically only a GUI showing Events, identified by a specific string
% describing the type of event observed within one specific ROI.It is not
% more than a helper function, collecting specific Data from a rather large
% Data set. Within this functions not written by me are used, published on
% mathworks.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Time_Analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @Time_Analyzer_OutputFcn, ...
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


% --- Executes just before Time_Analyzer is made visible.
function Time_Analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for Time_Analyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Time_Analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Time_Analyzer_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;



function judgeStr_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of judgeStr as text
%        str2double(get(hObject,'String')) returns contents of judgeStr as a double


% --- Executes during object creation, after setting all properties.
function judgeStr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)

%loading all jedge strings for every experiment in every folder and
%subfolder
fn=uigetdir;
list=getfn(fn,'judge.mat');
handles.list=list';
list=list';
handles.res={};
%sorting the data concerning event judge_string using a specific string
%given by the user
handles.events={};
m=1;
for k=1:size(list,1)
    load(list{k});
    dmy=list{k};
    for i=1:size(jstr,1)
        if strcmp(jstr{i},handles.judgeStr.String)
            handles.events{m,1}=list{k};
            handles.events{m,2}=i;
            m=m+1;
        end
    end
end
assignin('base','jdg',handles.events)

handles.text2.String=[num2str(m-1),' events found']

guidata(hObject, handles);

function eventNo_Callback(hObject, eventdata, handles)
% Used to surf between the events
% Plotting is always done using the median as requested by the user
k=str2num(handles.eventNo.String);

pn=handles.events{k,1};
dmy=pn(1:end-9);
load([dmy,'results.mat']);

i=handles.events{k,2};
dmyl=results.Trace_l(:,i)-results.Trace_BG_l-(0.146*(results.Trace_r(:,i)-results.Trace_BG_r));

axes(handles.axes1);
%plot(dmyl,'b.');
%hold on
plot(movmedian(dmyl,20),'b.');
%hold off
drawnow

axes(handles.axes2);
%plot(results.Trace_r(:,i)-results.Trace_BG_r);
%hold on
plot(movmedian(results.Trace_r(:,i)-results.Trace_BG_r,20),'b.');
%hold off
drawnow

% Saving the results on a user specified path
pn=handles.savePath.String;
idx_top=[];
idx_bottom=[];
if exist([pn,'\',handles.judgeStr.String,'_results.mat'])==2  
    dmyres=load([pn,'\',handles.judgeStr.String,'_results.mat']);
    res=dmyres.res;   
    assignin('base','res',res);
    f=str2num(handles.eventNo.String);
    % setting bottom axes handles.axes2 and top axes handles.axes1
    no=str2num(handles.eventNo.String);
    if res{no,7}<res{no,9} %assuming first rows are bottom graphs, which never fails        
        for g=1:size(res{no,4},2) %getting events for bottom graphs            
            idx_bottom=[idx_bottom; res{no,4}{1,g}.Position];%collecting all data indices
        end
        for g=1:size(res{no,5},2) %getting events for top graph            
            idx_top=[idx_top; res{no,5}{1,g}.Position];%collecting all data indices
        end                 
    elseif isnan(res{no,9})
        if res{no,7}<10
            for g=1:size(res{no,4},2) %getting events for bottom graphs            
                idx_bottom=[idx_bottom; res{no,4}{1,g}.Position];%collecting all data indices
            end
        else
            for g=1:size(res{no,4},2) %getting events for bottom graphs            
                idx_top=[idx_top; res{no,4}{1,g}.Position];%collecting all data indices
            end
        end
    else        
        for g=1:size(res{no,5},2) %getting events for bottom graphs            
            idx_bottom=[idx_bottom; res{no,5}{1,g}.Position];%collecting all data indices
        end
        for g=1:size(res{no,4},2) %getting events for top graph
            idx_top=[idx_top; res{no,4}{1,g}.Position];%collecting all data indices
        end
    end
    
    % User defines which graph shows which signal
    if strcmp(res{no,6},'vesicle')
        handles.bottom.Value=1;
    elseif strcmp(res{no,6},'membrane')
        handles.bottom.Value=2;
    end
    
    if strcmp(res{no,8},'membrane')
        handles.up.Value=1;
    elseif strcmp(res{no,8},'vesicle')
        handles.up.Value=2;
    end
        
    axes(handles.axes2);
    hold on
    if ~isempty(idx_bottom)
        plot(idx_bottom(:,1),idx_bottom(:,2),'rx','MarkerSize',12,'Linewidth',2);
    end
    hold off
    %assignin('base','idx_top',idx_top)
    %size(idx_top)
    axes(handles.axes1);
    hold on
    if ~isempty(idx_top)
        plot(idx_top(:,1),idx_top(:,2),'rx','MarkerSize',12,'Linewidth',2);
    end
    hold off
   
end
drawnow
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function eventNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% Using the DataTips the user is allowed to specifiy as much positions in
% the graphs as he wants. The number depends on the type of event is solely
% userchoice
h=Time_Analyzer;
cursorobj = datacursormode(h);
cursorobj.SnapToDataVertex = 'on'; % Snap to our plotted data, on by default
while true
    w = waitforbuttonpress;
    % waitforbuttonpress returns 0 with click, 1 with key press
    % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
    cursorobj.Enable = 'on'; % Turn on the data cursor, hold alt to select multiple points
    key = get(gcf,'currentcharacter');
    switch key
        case 116
            break
        otherwise
    end
end
cursorobj.Enable = 'off';
%assignin('base','cursorobj',cursorobj);
mypoints = getCursorInfo(cursorobj);
%sorting the Data, by axes origin, needs to be done by position as tag gets
%lost due to misfunction in matlab
sCh=[];
for p=1:size(mypoints,2)
    a=get(mypoints(p).Target,'parent');
    sCh=[sCh;a.Position(2)];
end

chNo=unique(sCh);
m=1;
n=1;
ch1={};
ch2={};
for p=1:size(mypoints,2)
    if sCh(p)==chNo(1) %niedrige Position, untere Axes, muss nichr sein wenn nur ein Graph gesetzt wurde
        ch1{m}=mypoints(p);
        m=m+1;
    elseif sCh(p)==chNo(2)
        ch2{n}=mypoints(p);
        n=n+1;
    end
end

if size(chNo,1)==1
    chNo(2)=NaN;
end

k=str2num(handles.eventNo.String);

datacursormode(h,'off');

% Saving the data, in case of re-evaluation the user is asked if he wants
% to overwrite.

pn=handles.savePath.String;

if exist([pn,'\',handles.judgeStr.String,'_results.mat'])==2  
    dmyres=load([pn,'\',handles.judgeStr.String,'_results.mat']);
    res=dmyres.res;    
    if isempty(res{k})==1
        res{k,1}=handles.events(k,1);
        res(k,2)=handles.events(k,2);
        res{k,3}=handles.judgeStr.String;
        res{k,4}=ch1;
        res{k,5}=ch2;
        res{k,6}=handles.bottom.String{handles.bottom.Value};
        res{k,7}=chNo(1);
        res{k,8}=handles.up.String{handles.up.Value};
        res{k,9}=chNo(2);
    else
        answer = questdlg('Data exists, Overwrite?','Yes','No');
        switch answer
            case 'Yes'
                res{k,1}=handles.events(k,1);
                res(k,2)=handles.events(k,2);
                res{k,3}=handles.judgeStr.String;
                res{k,4}=ch1;
                res{k,5}=ch2; 
                res{k,6}=handles.bottom.String{handles.bottom.Value};
                res{k,7}=chNo(1);
                res{k,8}=handles.up.String{handles.up.Value};
                res{k,9}=chNo(2);                    
            case 'No'
        end        
    end    
    save([pn,'\',handles.judgeStr.String,'_results.mat'],'res');
else    
    res{k,1}=handles.events(k,1);
    res(k,2)=handles.events(k,2);    
    res{k,3}=handles.judgeStr.String;
    res{k,4}=ch1;
    res{k,5}=ch2; 
    res{k,6}=handles.bottom.String{handles.bottom.Value};
    res{k,7}=chNo(1);
    res{k,8}=handles.up.String{handles.up.Value};
    res{k,9}=chNo(2);
    save([pn,'\',handles.judgeStr.String,'_results.mat'],'res');
end
guidata(hObject, handles);
eventNo_Callback(hObject, eventdata, handles)



function savePath_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function savePath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in up.
function up_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function up_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bottom.
function bottom_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function bottom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextEvent.
% For ease of use, button to jump to the next event and cleaning of the
% necessary structs and cells
function nextEvent_Callback(hObject, eventdata, handles)
handles.eventNo.String=num2str(str2num(handles.eventNo.String)+1);
guidata(hObject, handles);
eventNo_Callback(hObject, eventdata, handles)


% --- Executes on button press in trash.
function trash_Callback(hObject, eventdata, handles)
answer = questdlg('Trash Data?','Yes','No');
switch answer
    case 'Yes'
        disp('yes')
        pn=handles.savePath.String;
        dmyres=load([pn,'\',handles.judgeStr.String,'_results.mat']);
        res=dmyres.res;  
        k=str2num(handles.eventNo.String);
        res{k,1}=handles.events(k,1);
        res(k,2)=handles.events(k,2); 
        res{k,3}='trash';
        res{k,4}=[];
        res{k,5}=[];
        res{k,6}=[];
        res{k,7}=[];
        res{k,8}=[];
        res{k,9}=[];
        save([pn,'\',handles.judgeStr.String,'_results.mat'],'res');
    case 'No'
end
