function varargout = quickPlotter(varargin)
% QUICKPLOTTER M-file for quickPlotter.fig
%      QUICKPLOTTER, by itself, creates a new QUICKPLOTTER or raises the existing
%      singleton*.
%
%      H = QUICKPLOTTER returns the handle to a new QUICKPLOTTER or the handle to
%      the existing singleton*.
%
%      QUICKPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICKPLOTTER.M with the given input arguments.
%
%      QUICKPLOTTER('Property','Value',...) creates a new QUICKPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quickPlotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quickPlotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quickPlotter

% Last Modified by GUIDE v2.5 08-Oct-2015 11:21:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quickPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @quickPlotter_OutputFcn, ...
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


% --- Executes just before quickPlotter is made visible.
function quickPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quickPlotter (see VARARGIN)

% Choose default command line output for quickPlotter
global hFig;
global hAx;
handles.output = hObject;
hFig=[];
hAx=[];
% Update handles structure
guiData=handles;
guiData.hPlot=createNewFigure(hObject);
guiData.updateGui=@updateGui;
guiData = restoreGUI(guiData);

qpData.colNames=[];
qpData.rowNames=[];
qpData.x1Ind=1;
qpData.x2Ind=1;
qpData.yInd=1;
qpData.x1Function=[];
qpData.x2Function=[];
qpData.yFunction=[];
qpData.plottingComman=[];
qpData.lastHandleStatus.handle=[];
qpData.lastHandleStatus.value=[];
qpData.lastHandleStatus.string=[];
qpData.SelectionCommand=[];
% guiData.updateSPM=@updateSPM;

guiData.qpData=qpData;

set(guiData.PlottingTypesLstBx,'String', {'Line' 'Scatter' 'Histogram' 'Group Summary' 'Boxplot'});
set(guiData.PlottingOptionsLstBx,'String',{'Y' 'X1' 'X2' 'X1||Y' 'X2||Y' 'X1, Y' 'X2, Y' 'X1, X2' 'Command'});
set(guiData.FunctionsLstBx,'String',{'Command line'});
set(guiData.Function2LstBx,'String',{'Command line'});
set(guiData.CommandLstBx,'String',{'Command line'});

guidata(hObject, guiData);

% UIWAIT makes quickPlotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function guiData=restoreGUI(guiData)
guiData.statusMap=StatusMap('quickPlotter');
restoreFunction1Lst(guiData.statusMap);
restoreFunction2Lst(guiData.statusMap);
restoreCommandLst(guiData.statusMap);

function restoreFunction1Lst(map)
function restoreFunction2Lst(map)
function restoreCommandLst(map)

% --- Outputs from this function are returned to the command line.
function varargout = quickPlotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in RowNameLstBx.
function RowNameLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to RowNameLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RowNameLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RowNameLstBx
CommonMethods.getSelection(handles. RowNameLstBx, true);
guiData=guidata(hObject);
synchRowNumber(hObject, get(guiData.RowNameLstBx,'Value'));

% --- Executes during object creation, after setting all properties.
function RowNameLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RowNameLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in YLstBx.
function YLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to YLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns YLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from YLstBx
guiData=guidata(hObject);
synchRowNumber(hObject, get(guiData.YLstBx,'Value'));


% --- Executes during object creation, after setting all properties.
function YLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in X1LstBx.
function X1LstBx_Callback(hObject, eventdata, handles)
% hObject    handle to X1LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns X1LstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from X1LstBx
guiData=guidata(hObject);
synchRowNumber(hObject, get(guiData.X1LstBx,'Value'));


% --- Executes during object creation, after setting all properties.
function X1LstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X1LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in X2LstBx.
function X2LstBx_Callback(hObject, eventdata, handles)
% hObject    handle to X2LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns X2LstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from X2LstBx
guiData=guidata(hObject);
synchRowNumber(hObject, get(guiData.X2LstBx,'Value'));

% --- Executes during object creation, after setting all properties.
function X2LstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X2LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in X1ChoiceLstBx.
function X1ChoiceLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to X1ChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns X1ChoiceLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from X1ChoiceLstBx
CommonMethods.getSelection(handles. X1ChoiceLstBx, true);
updateX1(hObject);
updateX1List(hObject);
plotData(hObject);
% --- Executes during object creation, after setting all properties.
function X1ChoiceLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X1ChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in X2ChoiceLstBx.
function X2ChoiceLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to X2ChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns X2ChoiceLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from X2ChoiceLstBx
CommonMethods.getSelection(handles. X2ChoiceLstBx, true);
updateX2(hObject);
updateX2List(hObject);
plotData(hObject);

% --- Executes during object creation, after setting all properties.
function X2ChoiceLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X2ChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlottingOptionsLstBx.
function PlottingOptionsLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to PlottingOptionsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlottingOptionsLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlottingOptionsLstBx
guiData=guidata(hObject);
updatePlottingOptions(hObject);
subInds=get(guiData.SubplotLstBx,'Value');
len=length(subInds);
subs=guiData.subplots;
if(len>1)%multiple subplots are selected
    for i=1:len
        subs{subInds(i)}.updatePlottingOption(hObject);
    end
end
plotData(hObject);

% --- Executes during object creation, after setting all properties.
function PlottingOptionsLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlottingOptionsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateGui(hObject, data, colNames, rowNames, appendData)

if(~exist('appendData','var'))
    appendData=false;
end
guiData=guidata(hObject);
if(~exist('data','var'))
    if(isfield(guiData,'qpData'))
        qpdata=guiData.qpData;
        data=qpdata.data;
    else
        data=[];
    end
end
if(~exist('colNames','var'))
    if(isfield(guiData,'qpData'))
        qpdata=guiData.qpData;
        colNames=qpdata.colNames;
    else
        colNames=[];
    end
end
if(~exist('rowNames','var'))
    if(isfield(guiData,'qpData'))
        qpdata=guiData.qpData;
        rowNames=qpdata.rowNames;
    else
        rowNames=[];
    end
end


if(~isempty(data))
    dim=size(data);
    if(length(rowNames)~=dim(1))
        rowNames=arrayfun(@(x) num2str(x), 1:dim(1), 'UniformOutput', false);
    end
    if(length(colNames)~=dim(2))
        colNames=arrayfun(@(x) num2str(x), 1:dim(2), 'UniformOutput', false);
    end
end

if(appendData)
    data0=guiData.qpData.data;
    if(size(data0,1)==size(data,1))
        data=[data0 data];
        colNames=horzcat(guiData.qpData.colNames, colNames);
    end
end

if(isfield(guiData,'qpData'))
    if(isfield(guiData.qpData, 'data'))
        if(size(guiData.qpData.data,1)~=size(data, 1))
            set(guiData.SubplotLstBx, 'Value', 1);
        end
    end
end

guiData.qpData.data=data;
guiData.qpData.colNames=colNames;
guiData.qpData.rowNames=rowNames;
guidata(hObject, guiData);

updateFuncValues(hObject);

updateYChoiceLst(hObject,colNames);
updateX1ChoiceLst(hObject,colNames);
updateX2ChoiceLst(hObject,colNames);
updateRowName(hObject,rowNames);
set(guiData.X1LstBx,'Value', 1);
set(guiData.X2LstBx,'Value', 1);
set(guiData.YLstBx,'Value', 1);
updateX1(hObject);
updateX2(hObject);
updateY(hObject);
updateYList(hObject);
updateX1List(hObject);
updateX2List(hObject);
updatePlottingOptions(hObject);
plotData(hObject);
updateSubplotLstBx(hObject);

function updatePlottingOptions(hObject)
guiData=guidata(hObject);
id=get(guiData.PlottingOptionsLstBx, 'value');
options=get(guiData.PlottingOptionsLstBx, 'String');
guiData.qpData.plottingOption=options(id);
id=get(guiData.PlottingTypesLstBx, 'value');
options=get(guiData.PlottingTypesLstBx, 'String');
guiData.qpData.plottingType=options(id);
guiData.qpData.SelectionCommand=get(guiData.SelectionCommandTF, 'String');
guidata(hObject, guiData);

function updateRowName(hObject, names)
guiData=guidata(hObject);
set(guiData.RowNameLstBx, 'String', names);

function updateYChoiceLst(hObject,choices)
%choices=horzcat(choices, {'Function2'});
guiData=guidata(hObject);
oldChoices=get(guiData.YChoiceLstBx, 'String');
len=length(choices);
if(~sameChoices(choices, oldChoices))
    set(guiData.YChoiceLstBx,'String', choices);
    if(get(guiData.YChoiceLstBx,'Value')>len)
        set(guiData.YChoiceLstBx,'Value', 1);
    end
    set(guiData.YChoiceLstBx,'Max',len);
end

function updateX1ChoiceLst(hObject,choices)
guiData=guidata(hObject);
oldChoices=get(guiData.X1ChoiceLstBx, 'String');
len=length(choices);
if(~sameChoices(choices, oldChoices))
    set(guiData.X1ChoiceLstBx,'String', choices);
    if(get(guiData.X1ChoiceLstBx,'Value')>len)
        set(guiData.X1ChoiceLstBx,'Value', 1);
    end
    set(guiData.X1ChoiceLstBx,'Max',len);
end

function updateX2ChoiceLst(hObject,choices)
%choices=horzcat(choices, {'Function1'});
guiData=guidata(hObject);
oldChoices=get(guiData.X2ChoiceLstBx, 'String');
len=length(choices);
if(~sameChoices(choices, oldChoices))
    set(guiData.X2ChoiceLstBx,'String', choices);
    if(get(guiData.X2ChoiceLstBx,'Value')>len)
        set(guiData.X2ChoiceLstBx,'Value', 1);
    end
    set(guiData.X2ChoiceLstBx,'Max',len);
end

function IS=sameChoices(c1, c2)
if(~iscell(c1))
    c1={c1};
end
if(~iscell(c2))
    c2={c2};
end
len=length(c1);
if(len~=length(c2))
    IS=false;
    return;
end
for i=1:len
    if(~strcmp(c1{i}, c2{i}))
        IS=false;
        return;
    end
end
    
IS=true;

function updateX1List(hObject)
guiData=guidata(hObject);
X1=quickPlotterGuiDataHandler.getX1_qpData(guiData.qpData);
X1=arrayfun(@(x) num2str(x), X1, 'UniformOutput', false);

set(guiData.X1LstBx,'String', X1);


function updateX2List(hObject)
guiData=guidata(hObject);
X2=quickPlotterGuiDataHandler.getX2_qpData(guiData.qpData);
X2=arrayfun(@(x) num2str(x), X2, 'UniformOutput', false);

set(guiData.X2LstBx,'String', X2);

function updateYList(hObject)
guiData=guidata(hObject);
Y=quickPlotterGuiDataHandler.getY_qpData(guiData.qpData);
if(~isempty(Y))
    YL=arrayfun(@(x) num2str(x), Y, 'UniformOutput', false);
    guiData=guidata(hObject);
    set(guiData.YLstBx,'String',YL);
end

% qpData.colNames=[];
% qpData.rowNames=[];
% qpData.x1Ind=1;
% qpData.x2Ind=1;
% qpData.yInd=1;
% qpData.x1Function=[];
% qpData.x2Function=[];
% qpData.yFunction=[];
% qpData.plottingComman=[];
% qpData.lastHandleStatus.handle=[];
% qpData.lastHandleStatus.value=[];
% qpData.lastHandleStatus.string=[];
% qpData.selectionCommand=[];


function updateX1(hObject)
if(exist('hObject','var'))
    guiData=guidata(hObject);
    id=get(guiData.X1ChoiceLstBx,'Value');
    choices=get(guiData.X1ChoiceLstBx, 'String');
    if(max(id)>length(choices))
        id=1;
    end
    choice=choices{id};
    guiData.qpData.x1Ind=id;
    if(strcmp(choice, 'Function1'))
        functionChoices=get(guiData.Function1LstBx, String);
        id=get(guiData.Function1LstBx, String);        
        guiData.qpData.x1Function=functionChoices{id};
    else
        guiData.qpData.x1Function=[];
    end
    guidata(hObject,guiData);
end

function updateX2(hObject)
if(exist('hObject','var'))
    guiData=guidata(hObject);
    id=get(guiData.X2ChoiceLstBx,'Value');
    choices=get(guiData.X2ChoiceLstBx, 'String');
    if(max(id)>length(choices))
        id=1;
    end
    choice=choices{id};
    guiData.qpData.x2Ind=id;
    if(strcmp(choice, 'Function1'))
        functionChoices=get(guiData.Function1LstBx, String);
        id=get(guiData.Function1LstBx, String);        
        guiData.qpData.x2Function=functionChoices{id};
    else
        guiData.qpData.x2Function=[];
    end
    guidata(hObject,guiData);
end

function updateY(hObject)
if(exist('hObject','var'))
    guiData=guidata(hObject);
    id=get(guiData.YChoiceLstBx,'Value');
    choices=get(guiData.YChoiceLstBx, 'String');
    if(max(id)>length(choices))
        id=1;
    end
    choice=choices{id};
    guiData.qpData.yInd=id;
    if(strcmp(choice, 'Function2'))
        functionChoices=get(guiData.Function2LstBx, String);
        id=get(guiData.Function2LstBx, String);        
        guiData.qpData.yFunction=functionChoices{id};
    else
        guiData.qpData.yFunction=[];
    end
    guidata(hObject,guiData);
end

function synchRowNumber(hObject, value)
guiData=guidata(hObject);
set(guiData.X1LstBx,'Value',value);
set(guiData.X2LstBx,'Value',value);
set(guiData.YLstBx,'Value',value);
set(guiData.RowNameLstBx,'Value',value);

function hFig = createNewFigure(hObject, cloneSub)
if(~exist('cloneSub','var'))
    cloneSub=false;
end
hFig=figure();
if(cloneSub)
    copySubplots(hObject,hFig);
else
    subplot(1,1,1);
    CommonMethods.arrangeSubplots(hFig, [1 1], @registerGCA);
end

guiDataF=guidata(hFig);
guiDataF.creatorApp='quickPlotter';
guiDataF.creatorHandle=hObject;
guidata(hFig, guiDataF);

guiData=guidata(hObject);
guiData.hPlot=hFig;
if(isfield(guiData, 'hAxis'))
%    guiData=rmfield(guiData, 'hAxis');
end
guidata(hObject,guiData);
updateSubplotLstBx(hObject);

function plotData(hObject)
global hPlot;
hPlot=hObject;
guiData = guidata(hObject);

choice=get(guiData.PlottingOptionsLstBx, 'Value');
options=get(guiData.PlottingOptionsLstBx, 'String');
option=options{choice};
if(~isfield(guiData,'hPlot'))
     createNewFigure(hObject);
     guiData=guidata(hObject);
end

if(~CommonMethods.isFigureHandle(guiData.hPlot))
     createNewFigure(hObject);
     guiData=guidata(hObject);
end

guidata(hObject, guiData);

spInd=get(guiData.SubplotLstBx, 'Value');
spHandlers=CommonMethods.getSubplotHandlers(guiData.hPlot);
sphs=spHandlers(spInd);
if(length(sphs)==1)
    sph=sphs{1};%update(obj,xData,yData,hFig,xInds,yInds,title,xLabel,yLabel,option, style)
    sph.updateQPGuiHandler(hObject);
    sph.plot();
else
    for i=1:length(sphs)
        sph=sphs{i};
        sph.updateData(hObject);
        sph.plot();
    end
end

%title(option,'fontsize',18);
ax=gca;
guiData.hAxis=ax;
guiData.subplots=CommonMethods.getSubplotHandlers(guiData.hPlot);

if(get(guiData.UnifyYScalesCB,'Value'))
    unifyYScales(hObject);
end

if(get(guiData.UnifyXScalesCB,'Value'))
    unifyXScales(hObject);
end

guidata(hObject, guiData);
updateSubplotLstBx(hObject);

% --- Executes on selection change in YChoiceLstBx.
function YChoiceLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to YChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns YChoiceLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from YChoiceLstBx
CommonMethods.getSelection(handles. YChoiceLstBx, true);
updateY(hObject);
updateYList(hObject);
plotData(hObject);

% --- Executes during object creation, after setting all properties.
function YChoiceLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YChoiceLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SubplotTF_Callback(hObject, eventdata, handles)
% hObject    handle to SubplotTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SubplotTF as text
%        str2double(get(hObject,'String')) returns contents of SubplotTF as a double
global hFig;
global hAx;
guiData=guidata(hObject);
if(~CommonMethods.isFigureHandle(hFig))
    hFig=guiData.hPlot;
end
dim=CommonMethods.str2nums_decimal(get(guiData.SubplotTF,'String'));
figure(hFig);
CommonMethods.arrangeSubplots(hFig,dim,@registerGCA);
updateSubplotLstBx(hObject);

% --- Executes during object creation, after setting all properties.
function SubplotTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubplotTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GetCurrentAxisBT.
function GetCurrentAxisBT_Callback(hObject, eventdata, handles)
% hObject    handle to GetCurrentAxisBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hFig;
global hAx;
guiData=guidata(hObject);
guiData.hPlot=hFig;
guiData.hAxis=hAx;
guiData.subplots=CommonMethods.getSubplotHandlers(hFig);
guidata(hObject, guiData);
updateSubplotLstBx(hObject);

function registerGCA(varargin)
global hFig;
global hAx;
global hPlot;
h=gcf;
if(strcmp(CommonMethods.getCreatorApp(h),  'quickPlotter'))
    hFig=gcf;
    hAx=gca;
    cp=get(hAx,'CurrentPoint');
    if(ishandle(hPlot))
        GetCurrentAxisBT_Callback(hPlot);
    end
    hc=CommonMethods.getCreatorHandle(h);
    if(ishandle(hc))
        
        updateCurrentPoint(hc, hAx, cp);
    end
end

function updateCurrentPoint(hObject, hAx, cp)
guiData=guidata(hObject);
hFig=guiData.hPlot;
xy=cp(1,1:2);
subs=CommonMethods.getSubplotHandlers(hFig);
for i=1:length(subs)
    sub=subs{i};
    sub.markCurrentPoint(cp);
end
sub=CommonMethods.getOverlappingSubplot(get(hAx,'position'),subs);
currRow=sub.getCurrentRow(cp);
if(~isempty(currRow))
    synchRowNumber(hObject, currRow);
end


% --- Executes on selection change in SubplotLstBx.
function SubplotLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to SubplotLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SubplotLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SubplotLstBx

guiData=guidata(hObject);
if(CommonMethods.isFigureHandle(guiData.hPlot))
    subs=CommonMethods.getSubplotHandlers(guiData.hPlot);
    ID=get(guiData.SubplotLstBx,'Value');
    if(length(ID)==1)
        sub=subs{ID};
        ax=sub.hAxs;
        if(ishandle(ax))
            axes(ax);
            registerGCA();
            GetCurrentAxisBT_Callback(hObject);
        end
    end
end
% --- Executes during object creation, after setting all properties.
function SubplotLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubplotLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unifyYScales(hObject)
guiData=guidata(hObject);
guiDataF=guidata(guiData.hPlot);
subs=guiDataF.subHandlers;
ids=get(guiData.SubplotLstBx, 'Value');
Axs=[];
for i=1:length(ids)
    Axs=[Axs subs{i}.hAxs];
end
CommonMethods.unifyYScales(Axs);

function unifyXScales(hObject)
guiData=guidata(hObject);
guiDataF=guidata(guiData.hPlot);
subs=guiDataF.subHandlers;
ids=get(guiData.SubplotLstBx, 'Value');
Axs=[];
for i=1:length(ids)
    Axs=[Axs subs{i}.hAxs];
end
CommonMethods.unifyXScales(Axs);

function updateSubplotLstBx(hObject)
guiData=guidata(hObject);
subs=CommonMethods.getSubplotHandlers(guiData.hPlot);
if(length(subs)>0)
    list=cellfun(@(x) x.getDescription(), subs, 'UniformOutput', false);
    Axes=cellfun(@(x) num2str(x.hAxs(1)), subs, 'UniformOutput', false);
    st0=get(guiData.SubplotLstBx, 'String');
    set(guiData.SubplotLstBx, 'String', list);
    
    if(~CommonMethods.strcmp_cell(st0, list))
        set(guiData.SubplotLstBx,'Value', 1);
    end
    
    set(guiData.SubplotLstBx, 'max', length(subs));
    sub=subs{1};
    set(guiData.SubplotTF,'String',[num2str(sub.rows) 'X' num2str(sub.cols)]);

    if(length(get(guiData.SubplotLstBx,'Value'))==1)
        id=1;
        if(isfield(guiData, 'hAxis'))
            hAx=guiData.hAxis;
            id=find(ismember(Axes,num2str(hAx)));
        end
        if(~isempty(id))
            set(guiData.SubplotLstBx, 'Value',id);
        else
            set(guiData.SubplotLstBx, 'Value',1);
        end
    end
    set(guiData.CurrentFigerT, 'String', ['FIGURE: ' num2str(guiData.hPlot)]);
end


% --- Executes on selection change in PlottingTypesLstBx.
function PlottingTypesLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to PlottingTypesLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlottingTypesLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlottingTypesLstBx
guiData=guidata(hObject);
updatePlottingOptions(hObject);
subInds=get(guiData.SubplotLstBx,'Value');
len=length(subInds);
if(len>1)%multiple subplot selection
    subs=guiData.subplots;
    for i=1:len
        subs{subInds(i)}.updatePlottingType(hObject);
    end
end
plotData(hObject);


% --- Executes during object creation, after setting all properties.
function PlottingTypesLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlottingTypesLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UnifyYScalesCB.
function UnifyYScalesCB_Callback(hObject, eventdata, handles)
% hObject    handle to UnifyYScalesCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UnifyYScalesCB
unifyYScales(hObject)

% --- Executes on button press in UnifyXScalesCB.
function UnifyXScalesCB_Callback(hObject, eventdata, handles)
% hObject    handle to UnifyXScalesCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UnifyXScalesCB
unifyXScales(hObject)

function copyFigure(hObject)
    createNewFigure(hObject,true);
    plotData(hObject);
    
function copySubplots(hObject,hFig)
    guiData=guidata(hObject);
    subs=CommonMethods.getSubplotHandlers(guiData.hPlot);
    CommonMethods.makeSubplots(hFig,subs);
    guiData.subplots=CommonMethods.getSubplotHandlers(hFig);
    
        


% --- Executes on button press in NewFigureBT.
function NewFigureBT_Callback(hObject, eventdata, handles)
% hObject    handle to NewFigureBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
copyFigure(hObject);



function SelectionCommandTF_Callback(hObject, eventdata, handles)
% hObject    handle to SelectionCommandTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SelectionCommandTF as text
%        str2double(get(hObject,'String')) returns contents of SelectionCommandTF as a double
plotData(hObject);

% --- Executes during object creation, after setting all properties.
function SelectionCommandTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectionCommandTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YSelectionTF_Callback(hObject, eventdata, handles)
% hObject    handle to YSelectionTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YSelectionTF as text
%        str2double(get(hObject,'String')) returns contents of YSelectionTF as a double
plotData(hObject);


% --- Executes during object creation, after setting all properties.
function YSelectionTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YSelectionTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FunctionsLstBx.
function FunctionsLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to FunctionsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FunctionsLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FunctionsLstBx
guiData=guidata(hObject);
ind=get(guiData.FunctionsLstBx,'Value');
choices=get(guiData.FunctionsLstBx,'String');
choice=choices{ind};
if strcmp(get(handles.figure1,'SelectionType'),'open')
%    importAssociatedObject(hObject, el, alias);
else
%    set(guiData.FileNameAssociatedObjectTF,'String',el);
end

if(strcmpi(choice, 'Command line'))
    func=CommonMethods.getBaseWorkspaceVariable('Func');
    if(~isempty(func))
        names=get(guiData.FunctionsLstBx, 'String');
        if(~any(ismember(names,func)))
            names{end+1}=func;
%            set(guiData.FunctionsLstBx, 'String', horzcat(choices, {func}));
            set(guiData.FunctionsLstBx, 'String', names);
        end
    end
    updateGui(hObject);
end

function updateFuncValues(hObject)
%this function will evaluate func and add it to the data of qpdata
%func is a string, and it will added to the cellarray "columnNames" of
%qpdata
guiData=guidata(hObject);
funcs=get(guiData.FunctionsLstBx, 'String');
for i=2:length(funcs)
    guiData.qpData=quickPlotterGuiDataHandler.updateFuncValue(guiData.qpData,funcs{i});
end
guidata(hObject, guiData);

% --- Executes during object creation, after setting all properties.
function FunctionsLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FunctionsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Function2LstBx.
function Function2LstBx_Callback(hObject, eventdata, handles)
% hObject    handle to Function2LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Function2LstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Function2LstBx


% --- Executes during object creation, after setting all properties.
function Function2LstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Function2LstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CommandLstBx.
function CommandLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to CommandLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CommandLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CommandLstBx


% --- Executes during object creation, after setting all properties.
function CommandLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CommandLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in VoxelHistogramBT.
function VoxelHistogramBT_Callback(hObject, eventdata, handles)
% hObject    handle to VoxelHistogramBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=get(handles.RowNameLstBx, 'String');
fname=fname{get(handles.RowNameLstBx,'Value')};
vols=spm_read_vols(spm_vol(fname));
idx=find(~isnan(vols));
y=vols(idx);
figure();
hist(y);
parts=strsplit('/',fname);
set(get(gca,'title'),'String', parts{6}, 'fontSize', 24);
set(gca,'fontSize', 24);


% --- Executes on button press in ExportDataBT.
function ExportDataBT_Callback(hObject, eventdata, handles)
% hObject    handle to ExportDataBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%fname=spm_select([0 Inf],'image','Select an image file to load the volume');
guiData=guidata(hObject);
qpdata=guiData.qpData;
vars=evalin('base','line');
assignin('base', 'qpdata', qpdata);
cnames=CommonMethods.strrep_cell(qpdata.colNames,',', '_');
name=cnames{1};
p1=strfind(name,'(');
p2=strfind(name,')');
if(~isempty(p1)&&~isempty(p2))
    temp=name(p1+1:p2-1);
    temp=strrep(temp,'-', 'n');
    temp=strsplit(' ',temp);
    temp=CommonMethods.strcell2line(temp,'_');
    fname=['QPdata_' temp '.csv'];
end

[fname,path]=uiputfile('*.csv', 'Saving data as a csv file', fname);
fid=fopen(fullfile(path,fname),'wt');
if(fid<1)
    return;
end
rnames=qpdata.rowNames;
data=qpdata.data;
line=CommonMethods.strcell2line(horzcat({'filename'}, cnames),', ');
fprintf(fid,'%s\n',line);
for l=1:size(data,1)
    line=[rnames{l} ', ' CommonMethods.vec2line(data(l,:),', ')];
    fprintf(fid,'%s\n',line);
end
fclose(fid);


% --- Executes on button press in QPData2BaseBT.
function QPData2BaseBT_Callback(hObject, eventdata, handles)
% hObject    handle to QPData2BaseBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guiData=guidata(hObject);
qpdata=guiData.qpData;
CommonMethods.assignToBaseWrokspace('qpdata', qpdata);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
qpdata=CommonMethods.getBaseWorkspaceVariable('qpdata');
if(~isempty(qpdata))
    guiData=guidata(hObject);
    guiData.qpData=qpdata;
    guidata(hObject, guiData);
end
