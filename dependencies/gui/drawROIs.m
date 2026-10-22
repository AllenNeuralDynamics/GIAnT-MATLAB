classdef drawROIs < handle
    properties
        hFigPanel; %image panel
        hButtonPanel; %button panel
        hF; %figure handle
        hGrid;
        hAx; %axes handle
        hButton1; hButton1c;
        hButton2;
        hButton3;
        hButton4;
        hButtonCh;
        hEditCh;

        hROIs;
        roiData; %ROI information for immediate export
        hIm; %image object

        chDim = 3; %images are in XYC order
        %timeDim = 4; %images are in XYCT order
        showChannel = 1;

        %butterworth filter parameters for highpass
        % b2;
        % a2;

        IM;
        IMc; %correlation image; normalized to median 0 standard deviation 1
        IMsk; %skewness image; normalized to median 0 standard deviation 1
        IMact; %the image used to display activity
        IMavg; %average image, normalized to [0 1]
        bleach;
        aData; %alignment data
        dsFac; %downsampling factor of the image file, used to properly set the framerate;
        numChannels;
        fn;
        dr = [];
        drsave;
        fnsave;
        fnStem;

        defaultChLabels = {'iGluSnFr'; 'jRGECO'};
        Clevel = 3; %adjustable threshold for display of correlation image, in standard deviations
    end

    methods
        function obj = drawROIs(IM, dr, fn, shapedata)
            %IM is XYC
            obj.numChannels = size(IM,3);
            obj.dr = dr;
            obj.fn = fn;

            %gui creation
            obj.hF = uifigure('name', ['ROI TRACING - ' obj.fn], 'KeyPressFcn', @obj.kpf, ...
                'CloseRequestFcn', @obj.close); %closing the window keeps ROIs, same as Done
            obj.hGrid = uigridlayout(obj.hF);
            obj.hGrid.RowHeight = {50, '1x'};
            obj.hGrid.ColumnWidth = {'1x'};

            obj.hButtonPanel = uigridlayout(obj.hGrid);
            obj.hButtonPanel.RowHeight = {'1x'};
            obj.hButtonPanel.ColumnWidth = {100,100,100,100,100,100, 100,100,100};
            obj.hButton1 = uibutton(obj.hButtonPanel, 'Text', 'Draw Polygon', 'ButtonPushedFcn', @obj.drawPolygon);
            obj.hButton1c = uibutton(obj.hButtonPanel, 'Text', 'Draw Circle', 'ButtonPushedFcn', @obj.drawCircle);
            obj.hButton2 = uibutton(obj.hButtonPanel, 'Text', 'Save ROIs', 'ButtonPushedFcn', @obj.saveROIs);
            obj.hButton3 = uibutton(obj.hButtonPanel, 'Text', 'Load ROIs', 'ButtonPushedFcn', @obj.loadROIs);
            obj.hButton4 = uibutton(obj.hButtonPanel, 'Text', 'Done', 'ButtonPushedFcn', @obj.close);
            % for chIx = 1:obj.numChannels
            %     obj.hButtonCh(chIx) = uibutton(obj.hButtonPanel, 'Text', ['Show Ch' num2str(chIx) ':'], 'ButtonPushedFcn', @(x1,x2)(obj.setCh(chIx)));
            %     obj.hEditCh(chIx) = uieditfield(obj.hButtonPanel, 'Value', obj.defaultChLabels{chIx});
            % end

            obj.hAx = uiaxes(obj.hGrid);

            IMdisp = annotationDisplayImage(IM);

            obj.hIm = imshow(IMdisp,'parent', obj.hAx);
            set(obj.hIm, 'HitTest', 'off', 'pickableparts', 'none', 'clipping', 'on')
            set(obj.hAx, 'hittest', 'off', 'PickableParts', 'all')

            %plot outlines
            if nargin>3 && ~isempty(shapedata)
            hold(obj.hAx, 'on')
            for k = 1:length(shapedata)
                boundary = shapedata{k};
                plot(obj.hAx, boundary(:,2), boundary(:,1), 'r')
            end
            end
        end

        function setCh(obj,chIx)
            obj.showChannel = chIx;
            updateDisplay(obj);
        end

        function updateDisplay(obj)
            set(obj.hIm, 'CData',cat(3, obj.IMact(:,:,obj.showChannel)./obj.Clevel, obj.IMavg(:,:,obj.showChannel), obj.IMavg(:,:,obj.showChannel)))
        end

        function close(obj, evnt,~)
            roiData = {};
            for rix = 1:length(obj.hROIs)
                try
                    S = rmfield(get(obj.hROIs(rix)), {'Parent', 'Children', 'ContextMenu'});
                    S.mask = createMask(obj.hROIs(rix));
                    roiData = [roiData {S}];
                catch
                    %this object didn't exist
                end
            end
            obj.roiData = roiData;
            delete(obj.hF); %not close(), which would re-enter CloseRequestFcn
        end


        function saveROIs (obj, arg1, arg2)
            if isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_ROIs.mat']);
            if isempty(obj.fnsave) || ~any(obj.fnsave)
                return
            end
            roiData = {};
            for rix = 1:length(obj.hROIs)
                try
                    S = rmfield(get(obj.hROIs(rix)), {'Parent', 'Children', 'ContextMenu'});
                    roiData = [roiData {S}];
                catch
                    %this object didn't exist
                end
            end
            obj.roiData = roiData;
            save([obj.drsave filesep obj.fnsave], 'roiData');

            %save figure
            disp('saving figure')
            exportgraphics(obj.hAx, [obj.drsave filesep obj.fnStem '_FIGURE.pdf'], 'ContentType', 'vector');
        end

        function loadROIsDirect (obj, arg1)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
%             [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            load(arg1, 'roiData');

            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            disp(length(roiData))
            for rix = 1:length(roiData)
                %Lucas added this
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('ellipse was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.circle'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('circle was incomplete, creating a dud circle to delete later')
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', [5, 5], 'Label', 'dud');
                        else
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Radius', roiData{rix}.Radius, 'Label', roiData{rix}.Label);
                        end
                    
                    case 'images.roi.polygon'
                        if isempty(roiData{rix}.Position)
                            disp(rix)
                            disp('Polygon was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Label', roiData{rix}.Label);
                        end

                    otherwise
                        error('bad ROI data type');
                end

                fnms = fieldnames(roiData{rix});
                for ix =1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end
            end
        end

        function loadROIs (obj, arg1, arg2)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            load([dr filesep fn], 'roiData');

            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            disp(length(roiData))
            for rix = 1:length(roiData)
                %Lucas added this
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('ellipse was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.circle'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('circle was incomplete, creating a dud circle to delete later')
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', [5, 5], 'Label', 'dud');
                        else
                            obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Radius', roiData{rix}.Radius, 'Label', roiData{rix}.Label);
                        end
                    
                    case 'images.roi.polygon'
                        if isempty(roiData{rix}.Position)
                            disp(rix)
                            disp('Polygon was incomplete, creating a dud ellipse to delete later')
                            obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');                        
                        else
                            obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Label', roiData{rix}.Label);
                        end

                    otherwise
                        error('bad ROI data type');
                end

                fnms = fieldnames(roiData{rix});
                for ix =1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end


                cm = get(obj.hROIs(rix), 'ContextMenu');
                uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(rix)))); %add a 'Set Label' context menu
            end
        end

        function drawPolygon(obj, arg1, arg2) 
            L = length(obj.hROIs)+1;
            if L==1
                obj.hROIs = drawpolygon(obj.hAx);
            else
                obj.hROIs(L) = drawpolygon(obj.hAx);
            end
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function drawCircle(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            if L==1
                obj.hROIs = drawellipse(obj.hAx);
            else
                obj.hROIs(L) = drawellipse(obj.hAx);
            end
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function setLabel(obj, hROI)
            label = inputdlg('Set a label');
            if ~isempty(label)
                set(hROI, 'Label', label{1});
            end
        end
        function kpf(obj, arg1, arg2)
            switch arg2.Key
                case 'a'
                    @obj.drawPolygon;
                case 'c'
                    @obj.drawCircle;
                case 'period'
                    obj.Clevel = obj.Clevel/1.5;
                    obj.updateDisplay
                case 'comma'
                    obj.Clevel = obj.Clevel.*1.5;
                    obj.updateDisplay;
            end
        end
    end
end