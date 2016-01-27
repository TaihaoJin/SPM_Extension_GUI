classdef FourD_Vol_Handler < handle
    %FOURD_VOL_HANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coor;
        fname;
        name;
        hPlotter;
        vols;
        mat;
        type;
        rowNames;
        additionalNames;
        additionalData;
        storedCoors;
        errMsg;
    end
    
    methods
        function obj = FourD_Vol_Handler(fname, hPlotter)
            obj.fname=fname;
            obj.name=fname;
            if(~exist('hPlotter', 'var'))
                hPlotter=[];
            end               
            obj.hPlotter=hPlotter;
            [~, name, ext]=fileparts(fname);
            obj.storedCoors={};
            
            if(strcmpi(name, 'SPM')&&strcmpi(ext,'.mat'))
                obj.type='SPM';
                SPM=CommonMethods.spm_loadMat(fname);
                obj.vols=SPM.vols;
                dim=size(obj.vols);
                clear obj.SPM.vol;
                obj.mat=SPM.mat;
                clear SPM.mat;
                scans=SPM.xY.P;
                snames=CommonMethods.spm_getShortScanNames(scans);
                obj.rowNames=snames;
                obj.additionalData=SPM.xX.X;
                dim=size(SPM.xX.name);
                if(dim(1)>dim(2))
                    obj.additionalNames=SPM.xX.name';
                else
                    obj.additionalNames=SPM.xX.name;
                end
            else
                obj.type='4D';
                try 
                    obj.V=spm_vol(fname);
                    obj.vols=spm_read_vols(obj.V);
                    obj.mat=obj.V.mat;
                    obj.rowNames=arrayfun(@(x) num2str(x), 1:size(obj.vols,4));
                    obj.additionalNames={};
                    obj.additionalData=[];
                catch em
                    obj.type='invalid';
                    obj.errMsg=em;
                end
            end            
        end
        
        function updatePlotter(obj,mni,appendData)
            if(~exist('appendData','var'))
                appendData=false;
            end
            
            if(exist('mni','var'))
                obj.coor=mni;
            else
                mni=obj.coor;
            end
            
            if(~CommonMethods.isFigureHandle(obj.hPlotter))
                obj.buildPlotter();
            end
            
            [y,vname]=getVoxelValue(obj,mni);
            datat=zeros(length(y),length(obj.storedCoors)+1);
            colNamest=cell(1,length(obj.storedCoors)+1);
            datat(:,1)=y;
            colNamest{1}=vname;
           
            for i=1:length(obj.storedCoors)
                [y,vname]=getVoxelValue(obj,obj.storedCoors{i});
                datat(:,i+1)=y;
                colNamest{i+1}=vname;
            end
            
            datat=[datat obj.additionalData];
            colNamest=horzcat(colNamest, obj.additionalNames);
            guiData = guidata(obj.hPlotter);
            guiData.updateGui(obj.hPlotter, datat, colNamest, obj.rowNames, appendData);
        end
        
        function plotData(obj, data, colNames, rowNames)
            guiData = guidata(obj.hPlotter);
            guiData.updateGui(obj.hPlotter, data, colNames,rowNames);
        end
        
        function [datat, colNamest, rowNames]= getData(obj,mni)            
            if(exist('mni','var'))
                obj.coor=mni;
            else
                mni=obj.coor;
            end
            
            if(~CommonMethods.isFigureHandle(obj.hPlotter))
                obj.buildPlotter();
            end
            
            [y,vname]=getVoxelValue(obj,mni);
            datat=zeros(length(y),length(obj.storedCoors)+1);
            colNamest=cell(1,length(obj.storedCoors)+1);
            datat(:,1)=y;
            colNamest{1}=vname;
           
            for i=1:length(obj.storedCoors)
                [y,vname]=getVoxelValue(obj,obj.storedCoors{i});
                datat(:,i+1)=y;
                colNamest{i+1}=vname;
            end
            
            datat=[datat obj.additionalData];
            colNamest=horzcat(colNamest, obj.additionalNames);
            rowNames=obj.rowNames;
        end
        
        function [y, name]=getVoxelValue(obj,mni,center)
            %mni should be a 3xcols maxtrix. each clomn is one set of
            %coordinates.
            if(size(mni,1) ~=3||iscell(mni))
                if(size(mni,2)==3)
                    mni=mni';
                else
                    error('mni should be a 3xcols maxtrix for the function ''getVoxelValue(...)''');
                end
            end
            
            if(iscell(mni))
                mni=cell2mat(mni);
            end
            
            if(size(mni,2)>1)
                %center is the coordinates used to label the data
                voiInd=CommonMethods.mni2ind_cols(mni,obj.mat);
                y=arrayfun(@(x)squeeze(obj.vols(voiInd(1,x),voiInd(2,x), voiInd(3,x), :)), 1:size(voiInd,2),'UniformOutput',false);
                y=squeeze(y);
                y=cell2mat(y);
                y=mean(y,2);
                if(~exist('center','var'))
                    mni=mean(mni,2);
                else
                    mni=center;
                end
                cn='ClstrMean\_';
            else
                voiInd=CommonMethods.mni2ind(mni,obj.mat);
                y=squeeze(obj.vols(voiInd(1),voiInd(2), voiInd(3), :));
                cn='SVxl\_';
            end
            name=['(' sprintf('%3.0f',mni(1)) ' ' sprintf('%3.0f',mni(2)) ' ' sprintf('%3.0f', mni(3)) ')\_' cn obj.name];
        end
        
        function buildPlotter(obj)
            obj.hPlotter=quickPlotter();
        end
        
        function setName(obj,name)
            obj.name=name;
        end
        
        function storeCoordinates(obj)
            coors0=obj.storedCoors;
            coors=cell(1,(length(coors0)+1));
            coors{1}=obj.coor;
            for i=1:length(coors0)
                coors{i+1}=coors0{i};
            end
            obj.storedCoors=coors;
        end
        
        function clearStoredCoordinates(obj)
            obj.storedCoors={};
        end
    end
end

