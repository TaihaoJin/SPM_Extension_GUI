classdef quickPlotterGuiDataHandler <handle
    %QUICKPLOTTERGUIDATAHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        qpData;
        xType;
    end
    
    methods
        function obj = quickPlotterGuiDataHandler()
            obj.qpData=[];
            obj.xType='row';
        end       
        function update(obj,hObject)
            guiData=guidata(hObject);
            obj.qpData=guiData.qpData;
        end
        
        function updateData(obj,hObject)
            %update the data only. It assumes that the quickPlotter has
            %identical data structure as this object was updated last time.
            guiData=guidata(hObject);
            obj.qpData.data=guiData.qpData.qpData;
        end
        
        function updatePlottingOption(obj,hObject)
            guiData=guidata(hObject);
            options=get(guiData.PlottingOptionsLstBx, 'String');
            obj.qpData.plottingOption=options{get(guiData.PlottingOptionsLstBx, 'Value')};
        end
        
        function updatePlottingType(obj,hObject)
            guiData=guidata(hObject);
            types=get(guiData.PlottingTypesLstBx, 'String');
            obj.qpData.PlottingType=types{get(guiData.PlottingTypesLstBx, 'Value')};
        end
        
        function plot(obj,sph)
            CommonMethods.purgeOverlappingAxes(sph.hFig, sph.position)
            figure(sph.hFig);            
            sph.hAxs=axes('position', sph.position, 'ButtonDownFcn', sph.ButtonDownFcn);
            type=obj.qpData.plottingType;
            option=obj.qpData.plottingOption;
            indselection=quickPlotterGuiDataHandler.getSelection(obj.qpData);
            rows=length(indselection);
            if(strcmpi(option, 'y'))
                [Y,YL]=obj.getY_qpData(obj.qpData);
                X=(1:rows)';
                X=repmat(X,1,size(Y,2));
                XL='Index';
                obj.xType='row';
            elseif(strcmpi(option, 'x1'))
                [Y,YL]=obj.getX1_qpData(obj.qpData);
                X=(1:rows)';
                X=repmat(X,1,size(Y,2));
                XL='Index';
               obj.xType='row';
            elseif(strcmpi(option, 'x2'))
                X=(1:rows)';
                [Y,YL]=obj.getX2_qpData(obj.qpData);
                X=repmat(X,1,size(Y,2));
                XL='Index';
                obj.xType='row';
            elseif(strcmpi(option, 'x1, y'))                
                [Y,YL]=obj.getY_qpData(obj.qpData);
                [X,XL]=obj.getX1_qpData(obj.qpData);
                obj.xType='X1';
            elseif(strcmpi(option, 'x2, y'))
                [Y,YL]=obj.getY_qpData(obj.qpData);
                [X,XL]=obj.getX2_qpData(obj.qpData);
                obj.xType='X2';
            elseif(strcmpi(option, 'x1, x2'))
                [Y,YL]=obj.getX1_qpData(obj.qpData);
                [X,XL]=obj.getX2_qpData(obj.qpData);
                obj.xType='X1';
            elseif(strcmp(option,'Y Hist_Group'))
                xinds=obj.X1Inds;
                xinds=xinds(1);
                X1=obj.data(:,xinds);
                
                Y=obj.data(:,obj.YInds);
                Y=Y(X1>0);
                X=(1:length(Y))';
                
                XL='Index';
                YL=obj.colNames(obj.YInds);
                
                type='Histogram';
                obj.xType='Y';
            elseif(strcmpi(option, 'x1||y'))                
                [Y,YL1]=obj.getX1_qpData(obj.qpData);
                [Y2,YL2]=obj.getY_qpData(obj.qpData);                
                YL=horzcat(YL1,YL2);
                
                X=(1:rows)';
                X=repmat(X,1,size(Y,2));
                XL='Index';
                
                X2=(1:rows)';
                X2=repmat(X2,1,size(Y2,2));
                
                type='Plotyy';
                obj.xType='row';
            elseif(strcmpi(option, 'x2||y'))                
                [Y,YL1]=obj.getX2_qpData(obj.qpData);
                [Y2,YL2]=obj.getY_qpData(obj.qpData);                
                YL=horzcat(YL1,YL2);
                
                X=(1:rows)';
                X=repmat(X,1,size(Y,2));
                XL='Index';
                
                X2=(1:rows)';
                X2=repmat(X2,1,size(Y2,2));
                
                type='Plotyy';
                obj.xType='row';
            end
            
 %           axes(sph.hAxs(1));
            
            if(~exist('X', 'var'))
                X=1:length(Y);
                XL='Index';
            end
            if(~exist('Y', 'var'))
                [Y,YL]=obj.getY_qpData(obj.qpData);
            end
            if(strcmpi(type, 'Line'))
                plot(X,Y);
            elseif(strcmpi(type, 'Scatter'))
                scatter(X,Y);
                %scatter does not work on multiple column maxtrix X, Y
            elseif(strcmpi(type, 'Histogram'))
                hist(Y);
                XL=YL;
                YL='Counts';
                if(strcmp(option,'Y Hist_Group'))
                    YL='Counts\_Group';
                end
                obj.xType='Hist';
            elseif(strcmpi(type, 'Plotyy'))
                [Ax,~,~]=plotyy(X,Y,X2,Y2);
                set(Ax(1), 'ButtonDownFcn', sph.ButtonDownFcn);
                set(Ax(2), 'ButtonDownFcn', sph.ButtonDownFcn);
                ylabel(Ax(2),YL{2},'fontsize', 14);
                YL=YL{1};
                sph.hAxs=Ax;
            elseif(strcmpi(type, 'Summary'))

                if(isempty(sph.hAxs)||~all(ishandle(sph.hAxs)))
                     sph.hAxs=axes('position', sph.position, 'ButtonDownFcn', sph.ButtonDownFcn);
                end
                CommonMethods.plotSummary(sph.hAxs(1),Y);
                XL='Index';
                YL='Mean Sem';
            elseif(strcmpi(type, 'Boxplot'))
                if(isempty(sph.hAxs)||~all(ishandle(sph.hAxs)))
                     sph.hAxs=axes('position', sph.position, 'ButtonDownFcn', sph.ButtonDownFcn);
                end
                yInds=obj.qpData.yInd;
                xInds=obj.qpData.x1Ind;
                data=obj.qpData.data;
                inds=quickPlotterGuiDataHandler.getSelection(obj.qpData);                
                data=data(inds,:);
                grouping=true;
                for i=1:length(xInds)
                    xInd=xInds(i);
                    for j=1:length(yInds)
                        if(xInds==yInds(j))
                            grouping=false;
                            break;
                        end
                    end
                    if(~grouping)
                        break;
                    end
                    uq=unique(data(:,xInd));
                    if(length(uq)==2)
                        if(uq(1)~=0)
                            grouping=false;
                            break;
                        end
                    else
                        grouping=false;
                        break;
                    end
                end
                
                if(grouping)
                    Y1=cell(1,length(xInds));
                    for i=1:length(xInds)
                        Y1{i}=Y(data(:,xInds(i))~=0)
                    end
                else
                    Y1=data(:,yInds);
                end
                CommonMethods.boxplot_Cellarray(sph.hAxs(1),Y1);
                if(grouping)
                    XL='Group';
                else
                    XL='';
                end
            elseif(strcmpi(type, 'Group Summary'))
                if(isempty(sph.hAxs)||~all(ishandle(sph.hAxs)))
                     sph.hAxs=axes('position', sph.position, 'ButtonDownFcn', sph.ButtonDownFcn);
                end
                yInds=obj.qpData.yInd;
                xInds=obj.qpData.x1Ind;
                data=obj.qpData.data;
                inds=quickPlotterGuiDataHandler.getSelection(obj.qpData);
                data=data(inds,:);

                grouping=true;
                for i=1:length(xInds)
                    xInd=xInds(i);
                    for j=1:length(yInds)
                        if(xInds==yInds(j))
                            grouping=false;
                            break;
                        end
                    end
                    if(~grouping)
                        break;
                    end
                    uq=unique(data(:,xInd));
                    if(length(uq)==2)
                        if(uq(1)~=0)
                            grouping=false;
                            break;
                        end
                    else
                        grouping=false;
                        break;
                    end
                end
                
                if(grouping)
                    Y1=cell(1,length(xInds));
                    for i=1:length(xInds)
                        Y1{i}=Y(data(:,xInds(i))~=0)
                    end
                else
                    Y1=data(:,yInds);
                end
                CommonMethods.plotSummary(sph.hAxs(1),Y1);
                if(grouping)
                    XL='Group';
                else
                    XL='';
                end
                YL='Mean Sem';
            end
            xlabel(XL,'fontsize',14);
            ylabel(YL,'fontsize',14); 
            set(gca, 'ButtonDownFcn', sph.ButtonDownFcn);
            if(~strcmpi(type, 'Plotyy'))
                sph.hAxs=gca;
            end
            set(sph.hAxs,'position',sph.position);
        end
        
%         function inds=getSelectionX(obj,X)
%             selec=obj.XSelectionStr;
%             if(strcmpi(selec,'all'))
%                 inds=ones(1,size(obj.data,1));
%                 return
%             end
%             [lb, st1, sym, st2,hb]=CommonMethods.parseInequality(selec);
%             code=CommonMethods.getInequalityCode(st1, st2);
%             inds=arrayfun(@(x) CommonMethods.validInequality(lb,hb,x,code), X)';
%         end
%         function inds=getSelectionY(obj)
%             selec=obj.YSelectionStr;
%             if(strcmpi(selec,'all'))
%                 inds=ones(1,size(obj.data,1));
%                 return
%             end
%             [lb, st1, sym, st2,hb]=CommonMethods.parseInequality(selec);
%             nums=CommonMethods.str2nums_decimal(sym);
%             col=nums(1);
%             Y=obj.data(:,col);
%             code=CommonMethods.getInequalityCode(st1, st2);
%             inds=arrayfun(@(x) CommonMethods.validInequality(lb,hb,x,code),Y);
%         end
        
        function row = getCurrentRow(obj, cp)
            row=[];
            if(strcmpi(obj.xType,'row'))
                row=round(cp(1,1));
            end
        end
    end
    
    methods (Static=true)
        function [X1, labels]=getX1_qpData(qpData, evalSelection)
            formula=qpData.x1Function;
            if(~isempty(formula))
                X1=eval(formula);
                labels=arrayfun(@(x)[formula num2str(x)], 1:size(X1,2), 'UniformOutput', false);
            else
                X1=qpData.data(:,qpData.x1Ind);
                labels=qpData.colNames(qpData.x1Ind);
            end
            if(~exist('evalSelection','var'))
                evalSelection=true;
            end
            if(evalSelection)
                inds=quickPlotterGuiDataHandler.getSelection(qpData);
                X1=X1(inds);
            end
        end
        function [X2, labels]=getX2_qpData(qpData, evalSelection)
            formula=qpData.x2Function;
            if(~isempty(formula))
                varList=CommonMethods.getFunctionVariables(Formula);
                if(any(find(ismember(varList, 'X1'))))
                    X1=eval(formula);
                end
                X2=eval(formula);
                labels=arrayfun(@(x)[formula num2str(x)], 1:size(X2,2), 'UniformOutput', false);
            else
                X2=qpData.data(:,qpData.x2Ind);
                labels=qpData.colNames(qpData.x2Ind);
            end
            if(~exist('evalSelection','var'))
                evalSelection=true;
            end
            if(evalSelection)
                inds=quickPlotterGuiDataHandler.getSelection(qpData);
                X2=X2(inds);
            end
        end
        function [Y, labels]=getY_qpData(qpData, evalSelection)
            formula=qpData.yFunction;
            if(~isempty(formula))
                varList=CommonMethods.getFunctionVariables(formula);
                if(any(find(ismember(varList, 'X1'))))
                    X1=eval(formula);
                end
                if(any(find(ismember(varList, 'X2'))))
                    X2=eval(formula);
                end
                Y=eval(formula);
                labels=arrayfun(@(x)[formula num2str(x)], 1:size(Y,2), 'UniformOutput', false);
            else
                Y=qpData.data(:,qpData.yInd);
                labels=qpData.colNames(qpData.yInd);
            end
            if(~exist('evalSelection','var'))
                evalSelection=true;
            end
            if(evalSelection)
                inds=quickPlotterGuiDataHandler.getSelection(qpData);
                Y=Y(inds);
            end
        end
        function inds=getSelection(qpData)
            selec=qpData.SelectionCommand;
            Y=quickPlotterGuiDataHandler.getY_qpData(qpData, false);
            y=Y;
            X1=quickPlotterGuiDataHandler.getX1_qpData(qpData, false);
            x1=X1;
            X2=quickPlotterGuiDataHandler.getX2_qpData(qpData, false);
            x2=X2;
            if(strcmpi(selec,'all')||isempty(selec))
                inds=(1:size(qpData.data,1));
                return;
            else
                inds=find(eval(selec));
            end
        end
        
        function val=evalFunc(qpData,func)
            X1=quickPlotterGuiDataHandler.getX1_qpData(qpData);
            X2=quickPlotterGuiDataHandler.getX2_qpData(qpData);
            Y=quickPlotterGuiDataHandler.getY_qpData(qpData);
            data=qpData.data;
            colNames=qpData.colNames;
            rowNames=qpData.rowNames;
            
            val=eval(func);
        end
        
        function qpdata=updateFuncValue(qpData, func)
            qpdata=qpData;
            try
                temp=quickPlotterGuiDataHandler.evalFunc(qpData, func);
            catch ME
            end
            if(~exist('ME', 'var'))
                data=qpData.data;
                dim=size(data);
                idx=find(ismember(qpData.colNames,func));
                if(isempty(idx))
                    qpData.colNames{end+1}=func;
                    data=[data nan(dim(1),1)];
                    ind=dim(2)+1;
                else
                    ind=idx(1);
                end
                len=length(temp);
                len=min(len,dim(1));
                data(1:len,ind)=temp;
                qpData.data=data;
                qpdata=qpData;
            end
        end
    end
    
end

