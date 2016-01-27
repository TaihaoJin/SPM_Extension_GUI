function varargout = spm_extension(varargin)
% SPM_EXTENSION M-file for spm_extension.fig
%      SPM_EXTENSION, by itself, creates a new SPM_EXTENSION or raises the existing
%      singleton*.
%
%      H = SPM_EXTENSION returns the handle to a new SPM_EXTENSION or the handle to
%      the existing singleton*.
%
%      SPM_EXTENSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_EXTENSION.M with the given input arguments.
%
%      SPM_EXTENSION('Property','Value',...) creates a new SPM_EXTENSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_extension_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_extension_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spm_extension

% Last Modified by GUIDE v2.5 07-Dec-2015 14:32:39

% Begin initialization code - DO NOT EDIT
addpath(genpath('/home/JinT/matlab_projects/fMRI'));

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_extension_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_extension_OutputFcn, ...
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


% --- Executes just before spm_extension is made visible.
function spm_extension_OpeningFcn(hObject, eventdata, handles, varargin)
global guiData;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_extension (see VARARGIN)

% Choose default command line output for spm_extension
handles.output = hObject;
try
    X = load('TDdatabase');
    handles.DB = X.DB;
    handles.wholeMaskMNIAll = X.wholeMaskMNIAll;
catch
    errordlg('I can''t find TDdatabase.mat','TDdatabase not found');
end
% Update handles structure
guiData.statusMap=StatusMap('spm_extension');
guiData=handles;
guiData.setImageList=@setImageList;
guiData.update=@updateGUI;
guiData.setOpenImages=@setOpenImages;
guiData.updatePars=@updatePars;
guiData.spmStruct=[];
guiData.imgFileNameExt=get(guiData.FileNameAssociatedObjectTF,'String');
guiData.updatedMat=false;
guiData.statusMap=StatusMap('spm_extension');
guiData.ACPCAlignSuffix='AC_PC_Aligned';
guiData.SelectedImages=[];
set(guiData.ImageOperationLstBx,'String',{'Subtract Image', 'Apply Mask', 'Masking off Image', 'View SPM Results'});
set(guiData.ACPCOperationsLstBx,'String',{'Reset View' 'Zoom 40mm' 'Reset Origin' 'Align To Axis2' 'Align To 2-3 Plane' 'Save Aligned Image'});
%set(guiData.ImageHandlingLstBx,'String',{'Append Image' 'Close Image' 'Save Image' 'Intensity Profile' 'Voxel Histogram'});
set(guiData.ImageHandlingLstBx,'String',{'Append Image' 'Load Anatomic Mask' 'Close Image' 'Save Image' 'Intensity Profile' 'Voxel Histogram'});
guiData.temporaryImageFiles={};
anatomicRegionNames=fieldnames(guiData.wholeMaskMNIAll);
set(guiData.AnatomicRegionLstBx, 'String', anatomicRegionNames);
guidata(hObject, guiData);
updateRecentAssociatedObjectsLstBx(hObject);


% UIWAIT makes spm_extension wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_extension_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function ImageFileLB_Callback(hObject, eventdata, handles)
% hObject    handle to ImageFileLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageFileLB as text
%        str2double(get(hObject,'String')) returns contents of ImageFileLB as a double


% --- Executes during object creation, after setting all properties.
function ImageFileLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageFileLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FileNameAssociatedObjectTF_Callback(hObject, eventdata, handles)
% hObject    handle to FileNameAssociatedObjectTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileNameAssociatedObjectTF as text
%        str2double(get(hObject,'String')) returns contents of FileNameAssociatedObjectTF as a double
guiData=guidata(hObject);
guiData.imgFileNameExt=get(guiData.FileNameAssociatedObjectTF,'String');
if(isfield(guiData, 'associatedObjects'))
    ind=get(guiData.AssociatedObjectsLstBx,'value');
    aos=guiData.associatedObjects;
    ao=aos{ind};
    name=get(guiData.FileNameAssociatedObjectTF, 'String');
    ao.setName(name);
    names=get(guiData.AssociatedObjectsLstBx,'String');
    names{ind}=name;
    set(guiData.AssociatedObjectsLstBx,'String',names);
    
    if(isfield(guiData,'statusMap'))
        fname=ao.fname;
        srh=getRecentAssociatedObjectListHandler(hObject);
        srh.store(fname, name);
        guiData.statusMap.putSRHandler([],srh);
        updateRecentAssociatedObjectsLstBx(hObject);
    end
end

    function srh = getRecentAssociatedObjectListHandler(hObject)
        guiData=guidata(hObject);        
        key='AssociatedObjects';
        srh=guiData.statusMap.getSRHandler(key);  
        if(isempty(srh))
            srh=StringRankingHandler(key, 10, true);
        end
    
    
    function updateRecentAssociatedObjectsLstBx(hObject)
        guiData=guidata(hObject);
        srh=getRecentAssociatedObjectListHandler(hObject);
        [~,alias]=srh.getElements();
        set(guiData.RecentAssociatedObjectsLstBx, 'String', alias);
    
    
% --- Executes during object creation, after setting all properties.
function FileNameAssociatedObjectTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileNameAssociatedObjectTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateCoordinates_nmi(hObject,coorin)
guiData=guidata(hObject);    
app=getInteractingApp(hObject);
if(strcmpi(app,'spm'))
    guiData.spmStruct.centre=coorin;
elseif(strcmpi(app,'xjView'))
    guiData.xjViewStruct.currentxyz=coorin;
end
spm_orthviews('reposition',coorin,'fromSPMExt')
guidata(hObject, guiData);
updateGUI(hObject);

function output = updateGUI(hObject,varargin)
%This function updates the states of gui
guiData=guidata(hObject);    
set(guiData.ListingOptionListBx,'String',{'Opened images' 'Images to open'});
app=getInteractingApp(hObject);

if(strcmpi(app,'spm'))
    if(length(varargin)>=1)
        spmStruct=varargin{1};
        guiData.spmStruct=spmStruct;
    else
        spmStruct=guiData.spmStruct;
    end
    setOpenImages(hObject, CommonMethods.getOpenImages(spmStruct));
    
    imgIndex=CommonMethods.getCurrentVolumeIndex(spmStruct);%the index of the file that is associated with the clicked section view.
    V=spmStruct.vols{imgIndex};
    if(~isempty(V))
        % updateImageMatrix(hObject);
        
        options=get(guiData.ListingOptionListBx, 'String');
        index1=get(guiData.ListingOptionListBx,'Value');
        option=options{index1};
        switch (option)
            case 'Images to open'
            case 'Opened images'
                set(guiData.ImageFilesLstB, 'Value', imgIndex);
        end
        
        if(~isfield(V,'vol'))
            V.vol=spm_read_vols(spm_vol(V.fname));
            V.MaxStr={'Local Maxima List not Built'};
            V.MinStr={'Local Minima List not Built'};
            spmStruct.vols{imgIndex}=V;
            guiData.spmStruct=spmStruct;
        end
        
        guiData.updatedMat=false;
        coor=spmStruct.centre;
    end
elseif(strcmpi(app,'xjView'))
    
    if(length(varargin)>=1)
        xjViewStruct=varargin{1};
        guiData.xjViewStruct=xjViewStruct;
    else
        xjViewStruct=guiData.xjViewStruct;
    end

    guidata(hObject,guiData);
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    V=xjViewStruct.vols{imgIndex};  
    spmStruct=xjViewStruct.st;
    coor=xjViewStruct.currentxyz;
    if(isfield(V, 'clusterHandler'))
        stats=V.clusterHandler.stats;
        set(guiData.ClusterLstBx, 'String', stats.descriptions);
        cind=V.clusterHandler.getClusterIndx(coor);
        set(guiData.ClusterLstBx, 'Value', cind);
    end
end

if(~isempty(V))
    imgFile=V.fname;
    set(guiData.ImageFileLB,'String', imgFile);
    
    
    % updateImageMatrix(hObject);
    
    guiData.updatedMat=false;
    inds=CommonMethods.mni2ind(coor,V.mat);
    if(get(guiData.FindPeakRB, 'Value'))
        inds=CommonMethods.findPeak(V.vol, inds, 1);
        coor=CommonMethods.ind2mni_spmMat(V.mat, inds);
        spmStruct.centre=coor;
        guiData.spmStruct=spmStruct;
        spm_orthviews('update_st', spmStruct);
        spm_orthviews('redraw');
    end
    mmLB=sprintf('%6.2f %6.2f %6.2f', coor(1), coor(2), coor(3));
    volIndexes=inds;
    indLB=[num2str(volIndexes(1)) ' ' num2str(volIndexes(2)) ' ' num2str(volIndexes(3))];
    
    valLB=num2str(V.vol(volIndexes(1),volIndexes(2), volIndexes(3)));
    
    set(guiData.Coor_mmTF,'String', mmLB);
    set(guiData.Coor_voxelTF,'String', indLB);
    set(guiData.VoxelValueTF,'String', valLB);
    
    set(guiData.LocalMaximaLstBx,'String', V.MaxStr);
    set(guiData.LocalMinimaLstBx,'String', V.MinStr);
    [onelinestructure, cellarraystructure] = cuixuFindTDstructure(coor');
    set(guiData.AtlasLstBx, 'String', cellarraystructure);
    hObject=guiData.figure1;
    guidata(hObject,guiData);
    updatePlotter(hObject, coor);
end
updataAssociatedObjectListBox(hObject);
output=guiData;

% --- Executes on button press in DisplayImageBT.
function DisplayImageBT_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayImageBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
options=get(guiData.ListingOptionListBx, 'String');
index=get(guiData.ListingOptionListBx,'Value');
option=options{index};
app=getInteractingApp(hObject);

switch (option)
    case 'Images to open'
        if(isfield(guiData, 'ImagesToOpen'))
            if(strcmpi(app,'spm'))
                guiData=guidata(hObject);
                img=get(guiData.DisplayImageBT,'Value');
                spm_image('Init',img);
            end
        end
end

% --- Executes on selection change in ImageFilesLstB.
function ImageFilesLstB_Callback(hObject, eventdata, handles)
% hObject    handle to ImageFilesLstB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageFilesLstB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageFilesLstB
guiData=guidata(hObject);
options=get(guiData.ListingOptionListBx, 'String');
index=get(guiData.ListingOptionListBx,'Value');
option=options{index};
app=getInteractingApp(hObject);

switch (option)
    case 'Images to open'
        if(isfield(guiData, 'ImagesToOpen'))
            if(strcmpi(app,'spm'))
                Imgs=guiData.ImagesToOpen;
                if(iscell(Imgs))
                    if(~isempty(Imgs))
                        index=get(guiData.ImageFilesLstB,'Value');
                        set(guiData.SelectedImageIndexTF,'String', num2str(index));
                        spm_image('Init',img);
                        guiData=guidata(hObject);
                        CommonMethods.spm_setView(guiData.spmStruct, [0 0 0]);
                        CommonMethods.spm_setZoom(7);
                    end
                end
            end
        end
    case 'Opened images'   
        set(handles.ImageFilesLstB, 'Max', length(get(handles.ImageFilesLstB,'String')));
        if(strcmpi(app,'spm'))
            ind = get(guiData.ImageFilesLstB, 'Value');
            [~, names]=getOpenedImages(hObject);
            snames=get(handles.SelectedImageLstBx,'String')';
            set(handles.SelectedImageLstBx,'String', horzcat(snames, names(ind)));
            CommonMethods.spm_setCurrentVolumeIndex(guiData.spmStruct, ind);
            if(length(ind)==1)
                updateGUI(hObject);
            end
        elseif(strcmpi(app,'xjView'))
            updateGUI(hObject);
        end
end


% --- Executes during object creation, after setting all properties.
function ImageFilesLstB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageFilesLstB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setImageList(hObject,imgs)
guiData=guidata(hObject);
guiData.ImagesToOpen=imgs;
set(guiData.ImageFilesLstB,'Value',1);
set(guiData.ImageFilesLstB,'String',imgs);
guidata(hObject, guiData);

function setOpenImages(hObject,imgs)
guiData=guidata(hObject);
if(isfield(guiData,'OpenedImages'))
    imgs0=guiData.OpenedImages;
else
    imgs0=[];
end

updated=false;

if(length(imgs0)~=length(imgs))
    updated=true;
elseif(any(arrayfun(@(x) ~strcmp(imgs{x}, imgs0{x}), 1:length(imgs))))
    updated=true;
end

if (updated)
    guiData.OpenedImages=imgs;
    set(guiData.ImageFilesLstB,'Value',1);
    set(guiData.ImageFilesLstB,'String',CommonMethods.getFileNames(imgs));
    guidata(hObject, guiData);
end

function SelectedImageIndexTF_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedImageIndexTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SelectedImageIndexTF as text
%        str2double(get(hObject,'String')) returns contents of SelectedImageIndexTF as a double

%this function is written for adjusting "AC-PC" line using spm
guiData=guidata(hObject);
options=get(guiData.ListingOptionListBx, 'String');
index=get(guiData.ListingOptionListBx,'Value');
option=options{index};
app=getInteractingApp(hObject);

switch (option)
    case 'Images to open'
        if(isfield(guiData, 'ImagesToOpen'))
            if(strcmpi(app,'spm'))
                guiData=guidata(hObject);
                index=str2num(get(guiData.SelectedImageIndexTF,'String'));
                selectImageByIndex(hObject,index);
            end
        end
end


function selectImageByIndex(hObject,index)
guiData=guidata(hObject);
set(guiData.ImageFilesLstB,'Value',index);
imgs=get(guiData.ImageFilesLstB,'String');
img=imgs{index};
spm_image('Init',img);


% --- Executes during object creation, after setting all properties.
function SelectedImageIndexTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedImageIndexTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ListingOptionListBx.
function ListingOptionListBx_Callback(hObject, eventdata, handles)
% hObject    handle to ListingOptionListBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListingOptionListBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListingOptionListBx
persistent ImageFilesLstBValues;
persistent previousValue;
guiData=guidata(hObject);
value=get(guiData.ListingOptionListBx, 'Value');
if(isempty(ImageFilesLstBValues))
    ImageFilesLstBValues=ones(1,2);
end
if(isempty(previousValue))
    previousValue=1;
end

if(value~=previousValue)
    options=get(guiData.ListingOptionListBx, 'String');
    option=options{value};
    switch (option)
        case 'Images to open'
            if(isfield(guiData,'ImagesToOpen'))
                Imgs=guiData.ImagesToOpen;
            else
                Imgs={};
            end
            if(iscell(Imgs))
%                if(~isempty(Imgs))
                    set(guiData.ImageFilesLstB,'Value',ImageFilesLstBValues(value));
                    set(guiData.ImageFilesLstB,'String',Imgs);
%                end
            end
        case 'Opened images'
            Imgs=guiData.OpenedImages;
            if(iscell(Imgs))
                if(~isempty(Imgs))
                    set(guiData.ImageFilesLstB,'Value',ImageFilesLstBValues(value));
                    set(guiData.ImageFilesLstB,'String',Imgs);
                end
            end            
    end
else
    ImageFilesLstBValues(value)=get(guiData.ImageFilesLstB,'Value');
end
previousValue=value;

% --- Executes during object creation, after setting all properties.
function ListingOptionListBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListingOptionListBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Coor_mmTF_Callback(hObject, eventdata, handles)
% hObject    handle to Coor_mmTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coor_mmTF as text
%        str2double(get(hObject,'String')) returns contents of Coor_mmTF as a double
guiData=guidata(hObject);
str=get(guiData.Coor_mmTF,'String');
coor=CommonMethods.str2nums(' ', str)';
%spm_orthviews('reposition', coor);
% if(isfield(guiData, 'spmStruct'))
%     guiData.spmStruct.centre=coor;
%     guidata(hObject, guiData);
%     spm_orthviews('update_st', guiData.spmStruct);
% end
updateCoordinates_nmi(hObject,coor(1:3));

    function idx = getCurrentVoxelIndexes(hObject)
        guiData=guidata(hObject);
        str=get(guiData.Coor_voxelTF,'String');
        idx=CommonMethods.str2nums(' ', str)';
        
        
% --- Executes during object creation, after setting all properties.
function Coor_mmTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Coor_mmTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Coor_voxelTF_Callback(hObject, eventdata, handles)
% hObject    handle to Coor_voxelTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coor_voxelTF as text
%        str2double(get(hObject,'String')) returns contents of Coor_voxelTF as a double
guiData=guidata(hObject);
options=get(guiData.ListingOptionListBx, 'String');
id=get(guiData.ListingOptionListBx, 'Value');
option=options{id};
if(strcmp(option, 'Opened images'))
    inds=CommonMethods.str2nums_decimal(get(guiData.Coor_voxelTF, 'String'));
    app=CommonMethods.getSourceApp();
    imgIndex=get(guiData.ImageFilesLstB,'Value');
    if(strcmpi(app, 'SPM'))
        spmStruct=guiData.spmStruct;
        mat=spmStruct.vols{imgIndex}.mat;
        dim=spmStruct.vols{imgIndex}.dim;
    elseif(strcmpi(app,'xjView'))
        xjViewStruct=guiData.xjViewStruct;
        mat=xjViewStruct.vols{imgIndex}.mat;
        dim=xjViewStruct.vols{imgIndex}.dim;
    end
    
    if(length(inds)==3)
        inds = CommonMethods.getValidMatrixIndexes(dim, inds)
        coor=mat*([inds 1])';
        %spm_orthviews('reposition', coor(1:3));
        updateCoordinates_nmi(hObject,coor(1:3));
    end
end
% --- Executes during object creation, after setting all properties.
function Coor_voxelTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Coor_voxelTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VoxelValueTF_Callback(hObject, eventdata, handles)
% hObject    handle to VoxelValueTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VoxelValueTF as text
%        str2double(get(hObject,'String')) returns contents of VoxelValueTF as a double


% --- Executes during object creation, after setting all properties.
function VoxelValueTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoxelValueTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LocalMaximaLstBx.
function LocalMaximaLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to LocalMaximaLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LocalMaximaLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LocalMaximaLstBx
guiData = guidata(hObject);
value=get(guiData.LocalMaximaLstBx, 'Value');

app=getInteractingApp(hObject);
if(strcmpi(app,'spm'))
    index=CommonMethods.getCurrentVolumeIndex(guiData.spmStruct);
    V=guiData.spmStruct.vols{index};
else
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    V=guiData.xjViewStruct.vols{imgIndex};
end
mat=V.mat;
Inds=V.MaxPos(value,:);
Inds=[Inds [1]];
coor=mat*Inds';
%spm_orthviews('reposition',coor(1:3));
updateCoordinates_nmi(hObject,coor(1:3));

% --- Executes during object creation, after setting all properties.
function LocalMaximaLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LocalMaximaLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LocalMinimaLstBx.
function LocalMinimaLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to LocalMinimaLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LocalMinimaLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LocalMinimaLstBx
guiData = guidata(hObject);
value=get(guiData.LocalMinimaLstBx, 'Value');
app=getInteractingApp(hObject);
if(strcmpi(app,'spm'))
    index=CommonMethods.getCurrentVolumeIndex(guiData.spmStruct);
    V=guiData.spmStruct.vols{index};
else
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    V=guiData.xjViewStruct.vols{imgIndex};
end
mat=V.mat;
Inds=V.MinPos(value,:);
Inds=[Inds [1]];
coor=mat*Inds';
%spm_orthviews('reposition',coor(1:3));
updateCoordinates_nmi(hObject,coor(1:3));

% --- Executes during object creation, after setting all properties.
function LocalMinimaLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LocalMinimaLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function buildLocalExtrema(hObject)
guiData=guidata(hObject);

app=getInteractingApp(hObject);

if(strcmpi(app,'spm'))
    index=CommonMethods.getCurrentVolumeIndex(guiData.spmStruct);
    V=guiData.spmStruct.vols{index};
else
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    V=guiData.xjViewStruct.vols{imgIndex};
end

if(~isfield(V,'vol'))
    V.vol=spm_read_vols(spm_vol(V.fname));
end
[path,name,ext]=fileparts(V.fname);

stmfile=fullfile(path,[name '_xtm' '.csv']);

    vol=V.vol;
    if(length(size(vol))>3)
        vol=V.vol(:,:,:,1);
    end
    [Maxima,MaxPos,Minima,MinPos]=CommonMethods.localEtrema_clusterCenters(vol);

if(~CommonMethods.ValidAuxiliaryFile(V.fname, stmfile)||get(guiData.RecomputeLocalExtremaRB,'Value'))
    cprintf('blue','building local extrema list ... ...');
    vol=V.vol;
    if(length(size(vol))>3)
        vol=V.vol(:,:,:,1);
    end
    
    [Maxima,MaxPos,Minima,MinPos]=CommonMethods.localEtrema_clusterCenters(vol);
    cprintf('blue','finished local extrema list\n');
    exportXTM(stmfile,Maxima, MaxPos, Minima, MinPos);
else
    [Maxima,MaxPos,Minima,MinPos]=importXTM(stmfile);
end

ln=length(Minima);
lx=length(Maxima);
V.Maxima=Maxima;
V.MaxPos=MaxPos;
V.Minima=Minima;
V.MinPos=MinPos;
V.MaxStr=cell(1,lx);
V.MinStr=cell(1,ln);
len=max(ln,lx);
for i=1:len
    if(i<=lx)
        x=sprintf('%5d: %4d %4d %4d    %12.4d', i, MaxPos(i,1), MaxPos(i,2), MaxPos(i,3), Maxima(i));
        V.MaxStr{i}=x;
    end
    if(i<=ln)
        n=sprintf('%5d: %4d %4d %4d    %12.4d', i, MinPos(i,1), MinPos(i,2), MinPos(i,3), Minima(i));
        V.MinStr{i}=n;
    end
end

if(strcmpi(app,'spm'))
    guiData.spmStruct.vols{index}=V;
    spm_orthviews('update_st', guiData.spmStruct);
    updateGUI(hObject, guiData.spmStruct);
else
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    guiData.xjViewStruct.vols{imgIndex}=V;
    spm_orthviews('update_st', guiData.xjViewStruct.st);
    updateGUI(hObject, guiData.xjViewStruct);
end
guidata(hObject, guiData);

function exportXTM(stmfile, Maxima, MaxPos, Minima, MinPos)
nl=length(Minima);
nx=length(Maxima);
len=max(nl,nx);
table=-1*ones(len,8);
columnNames={'MaxI' 'MaxJ' 'MaxK' 'MaxV' 'MinI' 'MinJ' 'MinK' 'MinV'};
rowNames={};
for i=1:len
    if(i<=nx)
        table(i,1:3)=MaxPos(i,1:3);
        table(i,4)=Maxima(i);
    end
    if(i<=nl)
        table(i,5:7)=MinPos(i,1:3);
        table(i,8)=Minima(i);
    end
end
CommonMethods.exportTable_CSV(stmfile,columnNames,rowNames,table)

function [Maxima, MaxPos, Minima, MinPos] = importXTM(stmfile)
data=Importdata(stmfile);
data=data.data;
ni=data(:,1);
ni=data(ni>0);
len=length(ni);
MaxPos=data(1:len,1:3);
Maxima=data(1:len,4);
ni=data(:,5);
ni=data(ni>0);
len=length(ni);
MinPos=data(1:len,5:7);
Minima=data(1:len,8);


% --- Executes on button press in BuildLocalExtremaBT.
function BuildLocalExtremaBT_Callback(hObject, eventdata, handles)
% hObject    handle to BuildLocalExtremaBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
buildLocalExtrema(hObject);

function [onelinestructure, cellarraystructure] = cuixuFindTDstructure(pos)
% function [onelinestructure, cellarraystructure] = cuixuFindTDstructure(pos, DB, TALMNI)
%
% this function converts Talairach Daemon coordinate or MNI coordinate to a description of
% brain structure.
%
%   pos: the coordinate (MNI or Talairach coordinate, defined by TALMNI) of the position in the brain, in mm
%   DB: the database structure (see below for detailed explanation). If you don't input this argument, there
%   should be a file called 'TDdatabase.mat' in the current directory
%   TALMNI: 1 if pos is talairach coordinate, 0 if pos is MNI coordinate. 1
%   by default.
%
%   onelinestructure: a one-line description of the returned brain
%   structure
%   cellarraystructure: a cell array of size 5, each cell contains a string
%   describing the structure of the brain in a certain level.
%
%   Example:
%       [onelinestructure, cellarraystructure] = cuixuFindTDstructure([-24 24 48], DB, 1)
%       [onelinestructure, cellarraystructure] = cuixuFindTDstructure([-24 24 48])
%       then
%       onelinestructure = // Left Cerebrum // Frontal Lobe (L) // Middle Frontal Gyrus (L) // Gray Matter (L) // Brodmann area 8 (L)
%       cellarraystructure = 'Left Cerebrum'    'Frontal Lobe (L)'    [1x24 char]    'Gray Matter (L)'    'Brodmann area 8 (L)'
%
% Xu Cui
% 2004-6-28
%

%----------------------------------------------------------------------------------
% DB strcture
%----------------------------------------------------------------------------------
%-------------------------------
% Grid specification parameters:
%-------------------------------
% minX              - min X (mm)
% maxX              - max X (mm)
% voxX              - voxel size (mm) in X direction
% minY              - min Y (mm)
% maxY              - max Y (mm)
% voxY              - voxel size (mm) in Y direction
% minZ              - min Z (mm)
% maxZ              - max Z (mm)
% voxZ              - voxel size (mm) in Z direction
% nVoxX             - number of voxels in X direction
% nVoxY             - number of voxels in Y direction
% nVoxZ             - number of voxels in Z direction
%-------------------------------
% Classification parameters:
%-------------------------------
% numClass          - number of classification types
% cNames            - cNames{i}             - cell array of class names for i-th CT
% numClassSize      - numClassSize(i)       - number of classes for i-th CT
% indUnidentified   - indUnidentified(i)    - index of "indUnidentified" class for i-th CT
% volClass          - volClass{i}(j)        - number of voxels in class j for i-th CT
%
% data              - N x numClass matrix of referencies; let
%                       x y z coordinates in mm (on the grid) and
%                       nx = (x-minX)/voxX
%                       ny = (y-minY)/voxY
%                       nz = (z-minZ)/voxZ
%                       ind = nz*nVoxX*nVoxY + ny*nVoxX + nx + 1
%                       data(ind, i) - index of the class for i-th CT in cNames{i} to
%                                      which (x y z) belongs, i.e.
%                                      cNames{i}{data(ind, i)} name of class for i-th CT
%----------------------------------------------------------------------------------


%load('TDdatabase.mat');
%Modified by Taihao on 2/3/2015 to avoid possible confusion to other users.
%It was part of showsrs3
persistent wholeMaskMNIAll;
persistent DB;

if(isempty(DB))
    a=load('/home/kwj6/scripts/showsrs3/TDdatabase.mat');
    wholeMaskMNIAll=a.wholeMaskMNIAll;
    DB=a.DB;
end

if(size(pos,2)==4)
    pos=pos(1:3);
end
if(size(pos,1)>size(pos,2))
   pos=pos';
end

TALMNI = 0;

pos(:,1) = DB.voxX*round(pos(:,1)/DB.voxX);
pos(:,2) = DB.voxY*round(pos(:,2)/DB.voxY);
pos(:,3) = DB.voxZ*round(pos(:,3)/DB.voxZ);

min = [];
vox = [];
for(i=1:size(pos,1))
    min = [min; DB.minX DB.minY DB.minZ];
    vox = [vox; DB.voxX DB.voxY DB.voxZ];
end
n_pos = (pos - min)./vox;

nx = n_pos(:,1);
ny = n_pos(:,2);
nz = n_pos(:,3);
index = nz*DB.nVoxX*DB.nVoxY + ny*DB.nVoxX + nx + 1;
indMax = size(DB.data, 1);

onelinestructuretmp = [];
onelinestructure = [];
for(j=1:size(pos,1))
    onelinestructuretmp = [];
    for(i=1:DB.numClass)

        if (index(j) <= 0 | index(j) > indMax)
            ind(j) = DB.indUnidentified(i);
        else
            ind(j) = DB.data(index(j), i);
        end

        cellarraystructure{j,i} = DB.cNames{i}{ind(j)};
        onelinestructuretmp = [onelinestructuretmp ' // ' cellarraystructure{j,i}];
    end
    onelinestructure{j} = onelinestructuretmp;
end

return


% --- Executes on selection change in AtlasLstBx.
function AtlasLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to AtlasLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AtlasLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AtlasLstBx


% --- Executes during object creation, after setting all properties.
function AtlasLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AtlasLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindPeakRB.
function FindPeakRB_Callback(hObject, eventdata, handles)
% hObject    handle to FindPeakRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FindPeakRB


% --- Executes on button press in spmMatBT.
function spmMatBT_Callback(hObject, eventdata, handles)
% hObject    handle to spmMatBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
persistent spmViewer;

if(~CommonMethods.isFigureHandle(spmViewer))
    spmViewer=[];
end
if(isempty(spmViewer))
    spmViewer=spmMatViewer;
end

guiData=guidata(hObject);
st=guiData.spmStruct;
index=CommonMethods.getCurrentVolumeIndex(st)
vol=st.vols{index};
coor=st.centre;
ind=CommonMethods.mni2ind_spm(st, coor, index);
if(~isfield(vol,'SPM'))
    path=fileparts(vol.fname);
    matfile=fullfile(path,'SPM.mat');
    if(exist(matfile,'file'))
        a=load(matfile);
        SPM=a.SPM;
        SPM.filepath=path;
    else
        SPM=[];
    end
else
    SPM=vol.SPM;
end
if(~isempty(SPM))
    guiDatat=guidata(spmViewer);
    vol.SPM = guiDatat.updateSPM(spmViewer,SPM,ind);
    st.vols{index}=vol;
    guiData.spmStruct=st;
    guidata(hObject, guiData);
    spm_orthviews('update_st', st);
end


% --- Executes on selection change in AssociatedObjectsLstBx.
function AssociatedObjectsLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to AssociatedObjectsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AssociatedObjectsLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AssociatedObjectsLstBx
guiData=guidata(hObject);

if(isfield(guiData, 'associatedObjects'))
   ind=get(guiData.AssociatedObjectsLstBx, 'value');
   aoj=guiData.associatedObjects{ind};
   set(guiData.FileNameAssociatedObjectTF, 'string', aoj.fname);
end

% --- Executes during object creation, after setting all properties.
function AssociatedObjectsLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AssociatedObjectsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ClearAssociatedObjectBT.
function ClearAssociatedObjectBT_Callback(hObject, eventdata, handles)
% hObject    handle to ClearAssociatedObjectBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
if(isfield(guiData, 'associatedObjects'))
    aojs0=guiData.associatedObjects;
    ind=get(guiData.AssociatedObjectsLstBx, 'value');
    len=length(aojs0);
    if(len>1)
        aoj=aojs0{ind};
        objs=horzcat(aojs0{1:ind-1}, aojs0{ind+1:len});
        guiData.associatedObjects=objs;
        guidata(hObject, guiData);
        clear 'aoj';
    end
else
    set(guiData.AssociatedObjectsLstBx, 'String',{});
end
updataAssociatedObjectListBox(hObject);

    function updataAssociatedObjectListBox(hObject)
        guiData=guidata(hObject);
        names={};
        if(isfield(guiData, 'associatedObjects'))
            if(~isempty(guiData.associatedObjects))
                names=cellfun(@(x) x.name, guiData.associatedObjects, 'UniformOutput', false);
            end
        end
        set(guiData.AssociatedObjectsLstBx, 'String', names);
    
        
% --- Executes on button press in ImportAssociatedObjectBT.
function ImportAssociatedObjectBT_Callback(hObject, eventdata, handles)
% hObject    handle to ImportAssociatedObjectBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=spm_select([0 1],'Mat','Select an spm mat file or a 4D image file');
if(~isempty(fname))
    if(exist(fname,'file'))
        importAssociatedObject(hObject, fname);
    end
end

    function importAssociatedObject(hObject, fname, alias)
guiData=guidata(hObject);
if(isfield(guiData, 'associatedObjects'))
    associatedObjects=guiData.associatedObjects;
else
    associatedObjects={};
    set(guiData.AssociatedObjectsLstBx, 'String',{});
end

if(~exist('alias','var'))
    alias=fname;
end

set(guiData.FileNameAssociatedObjectTF, 'String', fname);
if(isfield(guiData, 'hPlotter'))
    hPlotter=guiData.hPlotter';
else
    hPlotter=[];
end

ID=-1;
for i=1:length(associatedObjects)
    fnamet=associatedObjects{i}.fname;
    if(strcmp(fname,fnamet))
        ID=i;
    end
end

if(ID>0)
    set(guiData.AssociatedObjectsLstBx, 'Value', ID);
else
    
    obj = FourD_Vol_Handler(fname, hPlotter);
    obj.name=alias;
    if(~strcmpi(obj.type, 'invalid'))
        %It is a new valid associated object
        
        srh=getRecentAssociatedObjectListHandler(hObject);
        srh.store(fname, alias);
        guiData.statusMap.putSRHandler([], srh);
        
        associatedObjects{end+1}=obj;
        set(guiData.AssociatedObjectsLstBx, 'Value', length(associatedObjects));
        guiData.associatedObjects=associatedObjects;
        guiData.hPlotter=obj.hPlotter;
        guidata(hObject,guiData);
    end
end
updateRecentAssociatedObjectsLstBx(hObject);
updataAssociatedObjectListBox(hObject);
size=length(get(guiData.AssociatedObjectsLstBx,'String'));
set(guiData.AssociatedObjectsLstBx,'Max',size);

function updatePlotter(hObject, mni)
guiData=guidata(hObject);
if(get(guiData.ClusterRB,'Value'))
    imgIndex=get(guiData.ImageFilesLstB, 'Value');
    V=guiData.xjViewStruct.vols{imgIndex};
    if(isfield(V,'clusterHandler'))
        handler=V.clusterHandler;
        [mnic, center, morphDist]=handler.getClusterMNIs(mni);
        if(~isempty(mnic))
            mni=mnic;
        end
    end
end
if(isfield(guiData, 'associatedObjects'))
    objs=guiData.associatedObjects;
    if(~isempty(objs))
        objs=objs(get(guiData.AssociatedObjectsLstBx, 'Value'));
        if(exist('mni', 'var'))
            obj0=objs{1};
            [data, colNames, rowNames]=obj0.getData(mni);
            
            for i=2:length(objs)
                obj=objs{i};
                [datat, colNamest, ~]=obj.getData(mni);
                data=[data datat];
                colNames=horzcat(colNames, colNamest);
            end
            obj0.plotData(data, colNames, rowNames);
        else
            objs{1}.updatePlotter();
        end
    end
end

function app = getInteractingApp(hObject)
%this function returns the label of the app that hObject is interacting
%with
guiData=guidata(hObject);
if(isfield(guiData,'interactingApp'))
    app=guiData.interactingApp;
else
    app=[];
end

function calling = callingSPM_orthviews(hObject, calling)
%this function set the status of hObject to indicate whether it has just
%called spm_orthviews.
guiData=guidata(hObject);
if(exist('calling','var'))
    guiData.callingSPM_orthview=calling;
    guidata(hObject,guiData);
else
    if(~isfield(guiData,'callingSPM_orthview'))
        calling=false;
    else
        calling=guiData.callingSPM_orthviews;
    end
end


% --- Executes on button press in StoreVoxelsBT.
function StoreVoxelsBT_Callback(hObject, eventdata, handles)
% hObject    handle to StoreVoxelsBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
if(isfield(guiData, 'associatedObjects'))
    for i=1:length(guiData.associatedObjects)
        guiData.associatedObjects{i}.storeCoordinates();
    end
    guidata(hObject, guiData);
    updatePlotter(hObject);
end

coorStrs=get(guiData.StoredCoordinatesLstBx,'String');
coorStrs{end+1}=get(guiData.Coor_mmTF,'String');
set(guiData.StoredCoordinatesLstBx,'String', coorStrs);
set(guiData.StoredCoordinatesLstBx,'Value', 1);

% --- Executes on button press in ClearVoxelsBT.
function ClearVoxelsBT_Callback(hObject, eventdata, handles)
% hObject    handle to ClearVoxelsBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
if(isfield(guiData, 'associatedObjects'))
    for i=1:length(guiData.associatedObjects)
        guiData.associatedObjects{i}.clearStoredCoordinates();
    end
    guidata(hObject, guiData);
    updatePlotter(hObject);
end
CommonMethods.removedSelection_LstBx(guiData.StoredCoordinatesLstBx);

    function updatePars(hObject, varargin)        
    %do nothing right now.
    
% --- Executes on button press in ClusterRB.
function ClusterRB_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ClusterRB


% --- Executes on selection change in ClusterLstBx.
function ClusterLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to ClusterLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ClusterLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClusterLstBx
guiData=guidata(hObject);    
set(guiData.ListingOptionListBx,'String',{'Opened images' 'Images to open'});
app=getInteractingApp(hObject);

if(strcmpi(app,'xjView'))
    imgIndex= get(guiData.ImageFilesLstB, 'Value');
    V=guiData.xjViewStruct.vols{imgIndex};    
    if(isfield(V, 'clusterHandler'))
        ind=get(guiData.ClusterLstBx, 'Value');
        stats=V.clusterHandler.stats;
        center=stats.centers(:,ind);
        mnis=stats.mnis{ind};
%         guiData.xjViewStruct.currentDisplayMNI=mat2cell(mnis);
        set(guiData.ClusterRB, 'Value', true);
        updateCoordinates_nmi(hObject, center);
        guidata(hObject, guiData);
    end
end



% --- Executes during object creation, after setting all properties.
function ClusterLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClusterLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
if(~isfield(guiData,'restoredGUI'))
    restoreGUI(hObject);
end

    function restoreGUI(hObject)
        guiData=guidata(hObject);
        statusMap=StatusMap('spm_extension');
        guiData.statusMap=statusMap;
        guiData.restoredGUI=true;
        guidata(hObject, guiData);
        

% --- Executes on selection change in RecentAssociatedObjectsLstBx.
function RecentAssociatedObjectsLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to RecentAssociatedObjectsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RecentAssociatedObjectsLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RecentAssociatedObjectsLstBx
guiData=guidata(hObject);
ind=get(guiData.RecentAssociatedObjectsLstBx,'Value');
strs=get(guiData.RecentAssociatedObjectsLstBx,'String');
alias=strs{ind};
srh=getRecentAssociatedObjectListHandler(hObject);
el=srh.getElement_alias(alias);
if strcmp(get(handles.figure1,'SelectionType'),'open')
    importAssociatedObject(hObject, el, alias);
else
    set(guiData.FileNameAssociatedObjectTF,'String',el);
end

% --- Executes during object creation, after setting all properties.
function RecentAssociatedObjectsLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RecentAssociatedObjectsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
saveGUI(hObject);
delete(hObject);
    function saveGUI(hObject)        
        guiData=guidata(hObject);
        if(isfield(guiData,'statusMap'))
            guiData.statusMap.save();
        end
        deletetemporaryFiles(hObject);
        
        


% --- Executes on button press in VoxelHistogamBT.
        function displayVoxelHistogam(hObject)
            % hObject    handle to VoxelHistogamBT (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            guiData=guidata(hObject);
            app=getInteractingApp(hObject);
            
            if(strcmpi(app,'spm'))
                spmStruct=guiData.spmStruct;
                setOpenImages(hObject, CommonMethods.getOpenImages(spmStruct));
                
                imgIndex=CommonMethods.getCurrentVolumeIndex(spmStruct);%the index of the file that is associated with the clicked section view.
                V=spmStruct.vols{imgIndex};
                if(~isempty(V))
                    % updateImageMatrix(hObject);
                    
                    options=get(guiData.ListingOptionListBx, 'String');
                    index1=get(guiData.ListingOptionListBx,'Value');
                    option=options{index1};
                    switch (option)
                        case 'Images to open'
                        case 'Opened images'
                            set(guiData.ImageFilesLstB, 'Value', imgIndex);
                    end
                    
                    if(~isfield(V,'vol'))
                        V.vol=spm_read_vols(spm_vol(V.fname));
                        V.MaxStr={'Local Maxima List not Built'};
                        V.MinStr={'Local Minima List not Built'};
                        spmStruct.vols{imgIndex}=V;
                        guiData.spmStruct=spmStruct;
                    end
                    
                    guiData.updatedMat=false;
                    coor=spmStruct.centre;
                end
            elseif(strcmpi(app,'xjView'))
                
                xjViewStruct=guiData.xjViewStruct;
                
                guidata(hObject,guiData);
                imgIndex= get(guiData.ImageFilesLstB, 'Value');
                V=xjViewStruct.vols{imgIndex};
                coor=xjViewStruct.currentxyz;
                if(isfield(V, 'clusterHandler'))
                    stats=V.clusterHandler.stats;
                    set(guiData.ClusterLstBx, 'String', stats.descriptions);
                    cind=V.clusterHandler.getClusterIndx(coor);
                    set(guiData.ClusterLstBx, 'Value', cind);
                end
            end
            
            vol=V.vol;
            if(V.dt(1)<=8)%interger voxel value
                idx=find(vol~=0);
            else
                idx=find(~isnan(vol));
            end
            y=vol(idx);
            figure();
            hist(y);
            stitle=CommonMethods.getFileName(V.fname);
            stitle=strrep(stitle,'_', '\_');
            title(stitle);
            set(gca, 'fontsize', 18);
            set(get(gca,'title'), 'fontsize', 24);
          
            


% --- Executes on button press in RecomputeLocalExtremaRB.
function RecomputeLocalExtremaRB_Callback(hObject, eventdata, handles)
% hObject    handle to RecomputeLocalExtremaRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RecomputeLocalExtremaRB

function [imgs, names]=getOpenedImages(hObject)
    imgs={};
    names={};
    guiData=guidata(hObject);
    app=getInteractingApp(hObject);
    
    if(strcmpi(app,'spm'))
        spmStruct=guiData.spmStruct;
        vols=spmStruct.vols;
    elseif(strcmpi(app,'xjView'))
        xjViewStruct=guiData.xjViewStruct;
        vols=xjViewStruct.vols;
    end
    
    for i=1:length(vols)
        vol=vols{i};
        if(~isempty(vol))
            imgs{end+1}=vol.fname;
            [~, name, ~]=fileparts(vol.fname);
            names{end+1}=[name];
        end
    end
        
function addImage(hObject, img)
    if(~exist('img','var'))
        [img, sts] = spm_select(1,'image','Select image');
    end
    imgs=getOpenedImages(hObject);
    imgs{end+1}=img;
    imgs=cellfun(@(x) spm_vol(x),imgs,'UniformOutput',false);
    imgs=cell2mat(imgs);
    spm_check_registration(imgs);
    addpath('/home/JinT/matlab_projects/fMRI/Scripts/Modified Scripts/spm');

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=get(handles.ImageFilesLstB, 'Value');
if(length(idx)~=2) return;end;

imgs=getOpenedImages(hObject);
imgs=imgs(idx);
V=spm_vol(imgs{1});
vols1=spm_read_vols(spm_vol(imgs{1}));
vols2=CommonMethods.getVoxelValues(imgs{1}, imgs{2});
vols2=vols2{1,1};
vold=vols1-vols2;
[path,name,ext]=fileparts(V.fname);
V.fname=fullfile(path,['temp' ext]);
spm_write_vol(V,vold);
addImage(hObject, V.fname);
delete(V.fname);

    function imgs = getSelectedImages(hObject)
        guiData=guidata(hObject);
        [oimgs, onames]=getOpenedImages(hObject);
        snames=get(guiData.SelectedImageLstBx,'String');
        if(~isempty(snames))
            imgs=cellfun(@(x) oimgs{ismember(onames,x)}, snames, 'UniformOutput', false);
        else
            imgs=[];
        end
        
% --- Executes on selection change in ImageOperationLstBx.
function ImageOperationLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to ImageOperationLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageOperationLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageOperationLstBx
options=get(handles.ImageOperationLstBx,'String');
operation=options{get(handles.ImageOperationLstBx,'Value')};

if(strcmp(operation, 'View SPM Results'))
    spm_results_ui();
end

imgs=getSelectedImages(hObject);
if(length(imgs)~=2) return;end;
V=spm_vol(imgs{1});
[path,name1,ext]=fileparts(V.fname);
[~,name2,~]=fileparts(imgs{2});
vols1=spm_read_vols(V);
vols2=CommonMethods.getVoxelValues(imgs{1}, imgs{2});
vols2=vols2{1,1};
fname='temp';
switch (operation)
    case 'Subtract Image'
        volsn=vols1-vols2;
        fname=[name1 '_subtract_' name2];
    case 'Apply Mask'
        if(V.dt<16)
            volsm=zeros(V.dim);
        else
            volsm=nan(V.dim);
        end
        volsn=vols1;
        m=vols2<0.5;
        volsn(m)=volsm(m);
        fname=[name1 '_maskOn_' name2];
    case 'Masking off Image'
        if(V.dt<16)
            volsm=zeros(V.dim);
        else
            volsm=nan(V.dim);
        end
        volsn=vols1;
        m=vols2>=0.5;
        volsn(m)=volsm(m);
        fname=[name1 '_maskOff_' name2];
    case 'View SPM Results'
        spm_results_ui();
end

V.fname=fullfile(path,[fname ext]);
spm_write_vol(V,volsn);
addImage(hObject, V.fname);
addTotemporaryFileList(hObject,V.fname);

    function addTotemporaryFileList(hObject,fname)
        guiData=guidata(hObject);
        guiData.temporaryImageFiles{end+1}=fname;
        guidata(hObject,guiData);
        
    function deletetemporaryFiles(hObject)
        guiData=guidata(hObject);
        fnames=guiData.temporaryImageFiles;
        cellfun(@(x) deleteImage(hObject,x), fnames);        

    function deleteImage(hObject,fname)
        closeImage(fname,hObject);
        delete(fname);

    function removeFromtemporaryFileList(hObject,fname)
        if(~iscell(fname))
            if(~ischar(fname))
                return;
            end
            fname={fname};
        end
        guiData=guidata(hObject);
        guiData.temporaryImageFiles=CommonMethods.removeStringEles(guiData.temporaryImageFiles, fname);
        guidata(hObject,guiData);
        
% --- Executes during object creation, after setting all properties.
function ImageOperationLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageOperationLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SelectedImageLstBx.
function SelectedImageLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedImageLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectedImageLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectedImageLstBx
idx=get(handles.SelectedImageLstBx,'Value');
strs=get(handles.SelectedImageLstBx,'String');
strs=CommonMethods.getStringCellarrayComparison(strs, strs(idx));
set(handles.SelectedImageLstBx,'Value', 1);
set(handles.SelectedImageLstBx,'String', strs{1});

% --- Executes during object creation, after setting all properties.
function SelectedImageLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedImageLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ACPCOperationsLstBx.
function ACPCOperationsLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to ACPCOperationsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ACPCOperationsLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ACPCOperationsLstBx
%'Reset View' 'Zoom 40mm' 'Reset Origin' 'Align To Axis2' 'Align To 2-3 Plane' 'Save Aligned Image'
options=get(handles.ACPCOperationsLstBx,'String');
operation=options{get(handles.ACPCOperationsLstBx,'Value')};

switch (operation)
    case 'Reset View'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            guiData=guidata(hObject);
            CommonMethods.spm_setView(guiData.spmStruct, [0 0 0]');
            CommonMethods.spm_setZoom(7);
        end
    case 'Zoom 40mm'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            CommonMethods.spm_setZoom(5);
        end
    case 'Reset Origin'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            guiData=guidata(hObject);
            st=CommonMethods.resetOrigin_spm(guiData.spmStruct);
            guiData.spmStruct=st;
            guiData.updatedMat=true;
            guidata(hObject,guiData);
            spm_orthviews('update_st', st);
            CommonMethods.spm_setView(st,[0 -23 0]);
            CommonMethods.spm_setZoom(5);
            %    spm_orthviews('Reposition',[0;0;0],'fromSPMExt');
        end
    case 'Align To Axis2'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            guiData=guidata(hObject);
            st=CommonMethods.toYAxis_spm(guiData.spmStruct);
            guiData.spmStruct=st;
            guiData.updatedMat=true;
            guidata(hObject,guiData);
            spm_orthviews('update_st', st);
            %    spm_orthviews('Reposition',[0;0;0],'fromSPMExt');
            CommonMethods.spm_setView(st,[0 20 0]);
            CommonMethods.spm_setZoom(7);
        end
    case 'Align To 2-3 Plane'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            guiData=guidata(hObject);
            st=CommonMethods.toYZPlane_spm(guiData.spmStruct);
            guiData.spmStruct=st;
            guiData.updatedMat=true;
            guidata(hObject,guiData);
            spm_orthviews('update_st', st);
            spm_orthviews('Reposition',[0;20;0],'fromSPMExt');
            %saving the AC_PC aligned image
            imgFile=get(guiData.ImageFileLB,'String');
            [path,name,ext]=fileparts(imgFile);
            name=[name '_' guiData.ACPCAlignSuffix];
            imgFile=fullfile(path,[name ext]);
            V=guiData.spmStruct.vols{1};
            V0=V;
            V.fname=imgFile;
            spm_write_vol(V,spm_read_vols(V0));
        end
    case 'Save Aligned Image'
        app=CommonMethods.getSourceApp();
        if(strcmpi(app,'SPM'))
            guiData=guidata(hObject);
            imgFile=get(guiData.ImageFileLB,'String');
            [path,name,ext]=fileparts(imgFile);
            name=[name '_' guiData.ACPCAlignSuffix];
            imgFile=fullfile(path,[name ext]);
            V=guiData.spmStruct.vols{1};
            V0=V;
            V.fname=imgFile;
            spm_write_vol(V,spm_read_vols(V0));
        end
end

% --- Executes during object creation, after setting all properties.
function ACPCOperationsLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACPCOperationsLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImageHandlingLstBx.
function ImageHandlingLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to ImageHandlingLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageHandlingLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageHandlingLstBx
%set(guiData.ImageHandlingLstBx,'String',{'Append Image' 'Close Image' 'Save Image' 'Intensity Profile' 'Voxel Histogram'});
options=get(handles.ImageHandlingLstBx,'String');
operation=options{get(handles.ImageHandlingLstBx,'Value')};
switch (operation)
    case 'Append Image'
        [img, sts] = spm_select(1,'image','Select image');
        imgs=getOpenedImages(hObject);
        imgs{end+1}=img;
        imgs=cellfun(@(x) spm_vol(x),imgs,'UniformOutput',false);
        imgs=cell2mat(imgs);
        spm_check_registration(imgs);
    case 'Load Anatomic Mask'
        imgs=getOpenedImages(hObject);
        img1=imgs{1};
        V=spm_vol(img1);
        
    case 'Close Image'
        imgs=getOpenedImages(hObject);
        idx=get(handles.ImageFilesLstB, 'Value');
        imgs=CommonMethods.removeCells(imgs, idx);
        imgs=cellfun(@(x) spm_vol(x),imgs,'UniformOutput',false);
        imgs=cell2mat(imgs);
        spm_check_registration(imgs);
    case 'Save Image'
        imgs=getOpenedImages(hObject);
        idx=get(handles.ImageFilesLstB, 'Value');
        arrayfun(@(x) removeFromtemporaryFileList(hObject,imgs{x}), idx);
    case 'Intensity Profile'
        imgs=getOpenedImages(hObject);
        idx=get(handles.ImageFilesLstB, 'Value');
        fname=imgs{idx(1)};
        V=spm_vol(fname);
        vols=spm_read_vols(spm_vol(fname));
        dim=size(vols);
        mat=V.mat;
        
        coori=mat*[0 0 0 1]';
        coorf=mat*[dim(1), dim(2), dim(3), 1]';
        coord=coorf-coori;
        
        idx = getCurrentVoxelIndexes(hObject);
        figure();
        
        [~, name, ~]=fileparts(fname);
        name=strrep(name,'_', '\_');
        subplot(2,2,1);
        len=dim(1);
        x=1:len;
        y=vols(:,idx(2),idx(3));
        x=x-idx(1);
        if(coord(1)<0) x=-1*x;end;
        plot(x,y);
        set(gca, 'fontsize', 14);
        title(name);
        xlabel('Left to Right (voxels)')
        ylabel('Intensity');
        
        subplot(2,2,2);
        len=dim(2);
        x=1:len;
        y=vols(idx(1),:,idx(3));
        x=x-idx(2);
        if(coord(2)<0) x=-1*x;end;
        plot(x,y);
        set(gca, 'fontsize', 14);
        title(name);
        xlabel('Posterior to Anterior (voxels)')
        ylabel('Intensity');        
        
        subplot(2,2,3);
        len=dim(3);
        x=1:len;
        y=squeeze(vols(idx(1),idx(2),:));
        x=x-idx(3);
        if(coord(3)<0) x=-1*x;end;
        plot(x,y);
        set(gca, 'fontsize', 14);
        title(name);
        xlabel('Inferior to Posterior (voxels)')
        ylabel('Intensity');
    case 'Voxel Histogram'
        displayVoxelHistogam(hObject);
end

    function closeImage(fname,hObject)
        imgs=getOpenedImages(hObject);
        idx=find(ismember(imgs,fname));
        if(isempty(idx))
            return;
        end
        imgs=CommonMethods.removeCells(imgs, idx);
        imgs=cellfun(@(x) spm_vol(x),imgs,'UniformOutput',false);
        imgs=cell2mat(imgs);
        spm_check_registration(imgs);

% --- Executes during object creation, after setting all properties.
function ImageHandlingLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageHandlingLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function SelectedImageLstBx_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to SelectedImageLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function MinIntensityTF_Callback(hObject, eventdata, handles)
% hObject    handle to MinIntensityTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinIntensityTF as text
%        str2double(get(hObject,'String')) returns contents of MinIntensityTF as a double


    function IntensityRangeSelection(hObject)
app=getInteractingApp(hObject);
guiData=guidata(hObject);
iI=str2num(get(guiData.MinIntensityTF, 'String'));
iF=str2num(get(guiData.MaxIntensityTF, 'String'));
if(strcmpi(app,'xjview'))
    if(isfield(guiData,'interactingAppHandle'))
        h=guiData.interactingAppHandle;
        guiData=guidata(h);
        guiData.CallBack_slider_TJ(guiData.slider,[],[iI iF]);
    end    
end

% --- Executes during object creation, after setting all properties.
function MinIntensityTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinIntensityTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MaxIntensityTF_Callback(hObject, eventdata, handles)
% hObject    handle to MaxIntensityTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxIntensityTF as text
%        str2double(get(hObject,'String')) returns contents of MaxIntensityTF as a double
IntensityRangeSelection(hObject);

% --- Executes during object creation, after setting all properties.
function MaxIntensityTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxIntensityTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StoredCoordinatesLstBx.
function StoredCoordinatesLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to StoredCoordinatesLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StoredCoordinatesLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StoredCoordinatesLstBx
guiData=guidata(hObject);
[idx, strs]=CommonMethods.getSelection(guiData.StoredCoordinatesLstBx);
set(guiData.Coor_mmTF,'String', strs{1});
Coor_mmTF_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function StoredCoordinatesLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StoredCoordinatesLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AnatomicMaskBT.
function AnatomicMaskBT_Callback(hObject, eventdata, handles)
% hObject    handle to AnatomicMaskBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function AnatomicRegionTF_Callback(hObject, eventdata, handles)
% hObject    handle to AnatomicRegionTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AnatomicRegionTF as text
%        str2double(get(hObject,'String')) returns contents of AnatomicRegionTF as a double
guiData=guidata(hObject);
set(guiData.AnatomicRegionLstBx, 'String', fieldnames(guiData.wholeMaskMNIAll));
strs=CommonMethods.searchStri_LstBx(guiData.AnatomicRegionLstBx, get(guiData.AnatomicRegionTF, 'String'));
set(guiData.AnatomicRegionLstBx, 'Value', length(strs));
set(guiData.AnatomicRegionLstBx, 'String', strs);

% --- Executes during object creation, after setting all properties.
function AnatomicRegionTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnatomicRegionTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in AnatomicRegionLstBx.
function AnatomicRegionLstBx_Callback(hObject, eventdata, handles)
% hObject    handle to AnatomicRegionLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnatomicRegionLstBx contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnatomicRegionLstBx
guiData=guidata(hObject);
[~,selection]=CommonMethods.getSelection(guiData.AnatomicRegionLstBx);
set(guiData.AnatomicRegionTF, 'String', selection{1});

% --- Executes during object creation, after setting all properties.
function AnatomicRegionLstBx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnatomicRegionLstBx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=spm_select();
fname=trim(fname(1,:));
[path, name, ext]=fileparts(fname);
fname=fullfile(path,'saved_coordinates.txt');
fid=fopen(fname,'wt');
guiData=guidata(hObject);
coors=get(guiData.StoredCoordinatesLstBx,'String');
for i=1:length(coors)
    fprintf(fid,'%s\n',coors{i});
end
fclose(fid);
