function varargout = dPETSTEP_4Dviewer(varargin)
    % GUI M-file for gui.fig
    %      GUI, by itself, creates a new GUI or raises the existings1
    %      singleton*.
    %
    %      H = GUI returns the handle to a new GUI or the handle to
    %      the existing singleton*.
    %
    %      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in GUI.M with the given input arguments.
    %
    %      GUI('Property','Value',...) creates a new GUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before gui_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to gui_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help gui

    % Last Modified by GUIDE v2.5 29-Jun-2017 09:39:33

    % Begin initialization code - DO NOT EDIT
    
gui_Singleton = 1;
    gui_State     = struct('gui_Name',   mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @gui_OpeningFcn, ...
                       'gui_OutputFcn',  @gui_OutputFcn, ...
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

end

% End initialization code - DO NOT EDIT

%==========================================================================
% --- Executes just before gui is made visible.
%==========================================================================
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to gui (see VARARGIN)
           
    % Check that there is some input
    if numel(varargin) == 0 
        disp('Invalid input to function!');
        delete(handles.figure1);
        return;
    end
       
    % GUI position on screen
    screenSize  = get(0, 'ScreenSize');
    set(hObject,'Units','pixels')
    guiSize     = get( hObject, 'OuterPosition');
    % Center window
%     set( hObject, 'OuterPosition', [ (screenSize(3)-guiSize(3))/2, (screenSize(4)-guiSize(4))/2, guiSize(3), guiSize(4) ] );
    set( hObject, 'OuterPosition', [0, 160, guiSize(3), guiSize(4) ] );
    set(hObject,'Units','pixels')

    % LOAD ICONS
    imageTrans = imread( fullfile(fileparts(mfilename('fullpath')),'image_transaxial.png') );  % Load image
    imageCor   = imread( fullfile(fileparts(mfilename('fullpath')),'image_coronal.png') );    % Load another image
    imageSag   = imread( fullfile(fileparts(mfilename('fullpath')),'image_sagittal.png') );
    imageTAC   = imread( fullfile(fileparts(mfilename('fullpath')),'image_TAC.png') );
    imageROIprof   = imread( fullfile(fileparts(mfilename('fullpath')),'image_ROIprofile.png') );
    imageROIHist   = imread( fullfile(fileparts(mfilename('fullpath')),'image_hist.png') );
    imshow((imageTrans),'Parent',handles.axesTrans);
    imshow((imageCor),'Parent',handles.axesCor);
    imshow((imageSag),'Parent',handles.axesSag);
    set( handles.TACTpushbutton,'CData',imageTAC);
    set( handles.ROIProfPushbutton,'CData',imageROIprof);
    set( handles.HistPushbutton,'CData',imageROIHist);
    
    % Current MATLAB version release date
    [~,handles.MATLABversion] = version; 
    
    % Other arguments than the image
    sizeImage = size(varargin{1}); if numel(sizeImage)<3; sizeImage(3)=1; end;
    if numel(varargin)>1
        for i=2:2:numel(varargin)-1
            if ischar(varargin{i})
                switch lower(varargin{i})
                    case {'frame'}
                        frame         = varargin{i+1};
                        frameLength   = diff(frame);
                        timeMidFrame  = diff(frame)/2 + frame(1:end-1);
                    case {'voxelsize'}
                        voxelSize     = varargin{i+1};
                    case {'sliceposition'}
                        slicePosition = varargin{i+1};
                    case {'dimunits'}
                        dimUnits = varargin{i+1};
                    case {'imageunit'}
                        imageUnit = varargin{i+1};
                end
            end
        end
    end 
    % Defaults if not set
    if ~exist('voxelSize','var')
        voxelSize = [1 1 1];
    end
    if ~exist('timeMidFrame','var')
        if numel(sizeImage)> 3 && sizeImage(4) > 1
            timeMidFrame = (0:1:sizeImage(4)-1)';
            frameLength  = ones(numel(timeMidFrame),1);
        else
            timeMidFrame = 0;
            frameLength  = 0;
        end
    end
    if ~exist('dimUnits','var')
        dimUnits = {'mm' 's'};
    end
    if ~exist('imageUnit','var')
        imageUnit = 'Bq/cc';
    end
    if ~exist('slicePosition','var')
        slicePosition = voxelSize(3)*(-(sizeImage(3)-1)/2:1:(sizeImage(3)-1)/2);
    end
    
    % Read data
    %----------------------------------------------------------------------
    handles.data            = varargin{1};
    handles.timeMidFrame    = timeMidFrame;
    handles.frameLength     = frameLength;
    handles.voxelSize       = voxelSize;
    handles.dimUnits        = dimUnits;
    handles.slicePosition   = slicePosition;
         
    % Clear some stuff
    clear slicePosition dimUnits voxelSize frameLength timeMidFrame
    
    % INITIALIZE 
    %----------------------------------------------------------------------

%     set(handles.figure1,'WindowScrollWheelFcn',@figScroll)
    set(handles.figure1,'WindowButtonDownFcn',@clickImageFunction)
    set(handles.figure1,'WindowButtonUpFcn',@unClickImageFunction)
    
    % Update handles structure
    guidata(hObject, handles);

    % ROIs
    handles.ROIs      = cell(0);
    handles.ROILegend = cell(0);
    handles.ROIColor  = {[0 0 1];[1 0 1];[0 1 0];[0 0.5 1];[1 1 0];[1 0 0];[0 1 1];...
                     [0.5 1 0];[1 0.5 0];[1 1 1];[0 1 0.5];[1 0 0.5];[0.5 0 1]};
%     handles.ROIColor  = {[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];...
%                      [0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]};
   
    % Projection view
    handles.proj     = 'tra';

    % Matrix and voxel size info
    set( handles.textVoxelSizeX,'String',[ num2str(handles.voxelSize(1)) ' ' handles.dimUnits{1} ]);
    set( handles.textVoxelSizeY,'String',[ num2str(handles.voxelSize(2)) ' ' handles.dimUnits{1} ]);
    set( handles.textVoxelSizeZ,'String',[ num2str(handles.voxelSize(3)) ' ' handles.dimUnits{1} ]);
    set( handles.textMatrixSizeX,'String',num2str(size(handles.data,1)));
    set( handles.textMatrixSizeY,'String',num2str(size(handles.data,2)));
    set( handles.textMatrixSizeZ,'String',num2str(size(handles.data,3)));
    
    handles.imageUnit = imageUnit;
    set( handles.textImageUnit,'String',handles.imageUnit );
    
    % Colormap
    set(handles.colormaps,'Value',1); % Jet
    colormaps_Callback(handles.colormaps, [], handles);
    handles.colorbar = colorbar(handles.axesColorbar);
    set(handles.colorbar, 'YAxisLocation','right');
    p = get(handles.axesColorbar,'Position');
    u = get(handles.axesColorbar,'Units');
    set(handles.colorbar,'Position',p,'Units',u);
    
    % SET ROI SLIDER MIN AND MAX VALUES
    %-----------------------------  
    handles.ROIwindow = 0.4;
    set(handles.ROIWindowSlider,'Min',0);
    set(handles.ROIWindowSlider,'Max',1);
    set(handles.ROIWindowSlider,'SliderStep',[0.05/1 0.5/1]);
    set(handles.ROIWindowSlider,'Value',handles.ROIwindow);
    
    % SET SLIDER STEP SIZE
    %---------------------
    % One click on arrow yields an increas of one slice, one click on the bar
    % yields an increase of 10 slices
    % If there is a third dimension...
    if size(handles.data,3)>1
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,3)-1)) 10/((size(handles.data,3)-1))]);
    %If there isn't...
    else
        set(handles.sliceSlider,'SliderStep',[1 1e4]);
    end
        
    % If there is a temporal dimension...
    if size(handles.data,4)>1
        set(handles.timeSlider,'SliderStep',[1/((size(handles.data,4)-1)) 10/((size(handles.data,4)-1))]);
    %If there isn't...
    else
        set(handles.timeSlider,'SliderStep',[1 1]);
    end

    % SET SLIDER MIN AND MAX VALUES
    %-----------------------------
    sliceSliderMin  = 1;
    sliceSliderMax  = size(handles.data,3);
    set(handles.sliceSlider,'Min',sliceSliderMin);
    set(handles.sliceSlider,'Max',sliceSliderMax);
    set(handles.textTotalSlice,'String',[ '/ ', num2str( size(handles.data,3) )] );

    timeSliderMin   = 1;
    timeSliderMax   = size(handles.data,4);
    set(handles.timeSlider,'Min',timeSliderMin);
    set(handles.timeSlider,'Max',timeSliderMax);
    set(handles.textTotalTime,'String',[ '/ ', num2str(timeSliderMax)] );

    % Inactivate sliders if max = min
    if sliceSliderMin == sliceSliderMax
        handles.sliceSlider.Enable = 'off'; 
%         handles.sliceSlider.Visible = 'off'; 
    end
    if timeSliderMin == timeSliderMax 
        handles.timeSlider.Enable = 'off';
%         handles.timeSlider.Visible = 'off'; 
    end

    lowSliderMin    = (min(handles.data(:)));
    lowSliderMax    = (max(handles.data(:)));
    if (lowSliderMax == lowSliderMin)
        step = get(handles.lowSlider,'sliderStep');
        lowSliderMax = lowSliderMin + step(1);
    end
    set(handles.lowSlider,'Min',lowSliderMin);
    set(handles.lowSlider,'Max',lowSliderMax);

    highSliderMin   = (min(handles.data(:)));
    highSliderMax   = (max(handles.data(:)));
    if (highSliderMax == highSliderMin)
        step = get(handles.highSlider,'sliderStep');
        highSliderMax = highSliderMin + step(1);
    end
    set(handles.highSlider,'Min',highSliderMin);
    set(handles.highSlider,'Max',highSliderMax);

    % Reset sliders upon open
    handles.sliderValueSpace = round( sliceSliderMax/2 );
    set(handles.sliceSlider,'Value' ,handles.sliderValueSpace);
    set(handles.sliceNumber,'String',num2str(handles.sliderValueSpace));
    sliceUnit = [ num2str( handles.slicePosition(handles.sliderValueSpace) ) ' ' handles.dimUnits{1} ];
    set(handles.sliceUnitText,'String', sliceUnit );
    
    handles.sliderValueTime = timeSliderMin;
    set(handles.timeSlider,'Value' ,timeSliderMin);
    set(handles.timeNumber,'String',num2str(timeSliderMin));
    timeUnit                = [ num2str( handles.timeMidFrame( handles.sliderValueTime ) ) ' ' handles.dimUnits{2} ];
    set(handles.timeUnitText,'String', timeUnit );
    
    % Min value of threshold slider
    handles.lowLevel        = lowSliderMin;
    set(handles.lowSlider,'Value' ,lowSliderMin);
    set(handles.lowNumber,'String',num2str(lowSliderMin));
    % Max value of threshold slider
    handles.highLevel       = highSliderMax;
    set(handles.highSlider,'Value' ,handles.highLevel);
    set(handles.highNumber,'String',num2str(handles.highLevel));
 
    % Sizes of the shown 2D data
    handles.xMin    = 1;
    handles.xMax    = size(handles.data,2);
    handles.yMin    = 1;
    handles.yMax    = size(handles.data,1);
    
    set( handles.textLineProfWidth,'String', ['1 / ' num2str(handles.xMax)] );

    % Max and min of image and slice, in information text box
    updateMinMaxText(handles);
    
    % Press "tra" button
    handles         = traPushbutton_Callback(handles.traPushbutton, [], handles);

    % Colorbar
    if (datenum(handles.MATLABversion)<datenum('October 3, 2014')) %Colorbar changed after R2014b (released 3/10-2014) 
        set( handles.colorbar, 'YLim',[handles.lowLevel handles.highLevel] );
        hCBImg = findobj(handles.colorbar,'type','image');
        set( hCBImg, 'YData',[handles.lowLevel handles.highLevel] );
    else
        set( handles.colorbar, 'Limits',[handles.lowLevel handles.highLevel] );
%         caxis([handles.lowLevel handles.highLevel])
        colormap( handles.colorbar, handles.colormaps.String{get(handles.colormaps,'Value')} );
        colormap( handles.axesColorbar, handles.colormaps.String{get(handles.colormaps,'Value')} );
        caxis(handles.axesColorbar,[handles.lowLevel handles.highLevel]);
    end
    
    % ROIs as input argument
    j=1;
    if numel(varargin)>1
        for i=2:2:numel(varargin)-1
            if ~ischar(varargin{i})
                sizeROI = size( varargin{i} );
                if numel(sizeROI)<=3 %2D or 3D ROI
                    
                    % Check that ROI have the right matrix size
                    if ~isequal([size(varargin{i},1),size(varargin{i},2),size(varargin{i},3)], [size(handles.data,1),size(handles.data,2),size(handles.data,3)])
                        fprintf('ROI "%s" has the wrong matrix dimension (%ix%ix%i). Should be (%ix%ix%i).\nExiting...\n\n',varargin{i+1},size(varargin{i},1),size(varargin{i},2),size(varargin{i},3),size(handles.data,1),size(handles.data,2),size(handles.data,3) );
                        delete(handles.figure1);
                        return;
                    end
                    
                    handles.ROIs{j} = varargin{i};

                    %Name
                    name = varargin{i+1};
                    set(handles.textROIName, 'String', name );
                    
                    handles = addROICont(handles);
                    
                    % Assume input ROI same orientation as input image data
                    if ~( size(handles.ROIs{j},1) == size(handles.data,1) &&...
                            size(handles.ROIs{j},2) == size(handles.data,2) &&...
                            size(handles.ROIs{j},3) == size(handles.data,3) )
                        handles.ROIs{j} = [];
                    end
                    j=j+1;
                elseif numel(sizeROI)>3 %ROI more than 3 dimensions not possible
                    fprintf('ROI "%s" has too many matrix dimensions. Should be 2 or 3.\nExiting...\n\n',varargin{i+1});
                    delete(handles.figure1);
                    return;
                end
            end
        end
    end

    % Choose default command line output for gui
    handles.output  = hObject;

    % Update handles structure
    guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
%     varargout{1} = handles.output;
end

% --- Executes on button press in traPushbutton.
function varargout = traPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to traPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    handles.proj = 'tra'; %Transaxial
    set(handles.projText,'String','Transaxial (YX slice)'); %'XY slice');

    %Delete line profile if any
    delete( findall( handles.axes1,'type','line') );
    set( handles.textLineProfWidth,'Enable', 'off' );
    set( handles.sliderLineProfWidth,'Enable', 'off' );
    set( handles.textLineProfText,'Enable', 'off' );
    set( handles.textLineProfWidth,'String', [num2str( get(handles.sliderLineProfWidth,'Value') )] );
    
    handles.sliderValueSpace     = round( get(handles.sliceSlider,'Value') );
    if handles.sliderValueSpace > size(handles.data,3)
        handles.sliderValueSpace = size(handles.data,3);
    end
          
    % Set slider values
    % If the data is 3D...
    if size(handles.data,3)>1
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,3)-1)) 10/((size(handles.data,3)-1))]);
    %If it isn't...
    else
        set(handles.sliceSlider,'SliderStep',[1 1e4]);
%         set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,3))) 10/((size(handles.data,3)))]);
    end
%     set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,3)-1)) 10/((size(handles.data,3)-1))])
    set(handles.sliceSlider,'Max',size(handles.data,3));
    set(handles.textTotalSlice,'String',[ '/ ', num2str( size(handles.data,3) )] );
    set(handles.sliceSlider,'value', handles.sliderValueSpace);
    sliceSlider_Callback( handles.sliceSlider,[],handles);
    
    % Inactivate sliders if max = min
    if size(handles.data,3)==1
        handles.sliceSlider.Enable = 'off';
    else
        handles.sliceSlider.Enable = 'on';
    end

    % Update text info at top
    updateMinMaxText(handles);
    
    handles.xMax = size(handles.data,2);
    handles.yMax = size(handles.data,1);

    % Choose axes to draw in
    axes(handles.axes1);
   
    % Plot image
    imshow(squeeze(handles.data(:,:,round(handles.sliderValueSpace),round(handles.sliderValueTime))),[handles.lowLevel handles.highLevel]);
    daspect( [ handles.voxelSize(1) handles.voxelSize(2) 1] );
    colormaps_Callback(handles.colormaps, [], handles)  
    if (datenum(handles.MATLABversion)<datenum('October 3, 2014')) %Colorbar changed after R2014b (released 3/10-2014) 
        set( handles.colorbar, 'YLim',[handles.lowLevel handles.highLevel] );
        hCBImg = findobj(handles.colorbar,'type','image');
        set( hCBImg, 'YData',[handles.lowLevel handles.highLevel] );
    else
        set( handles.colorbar, 'Limits',[handles.lowLevel handles.highLevel] );
        caxis(handles.axesColorbar,[handles.lowLevel handles.highLevel]);
    end
    axis off

    % Text in axes
    legendStr = get( handles.ROINumberMenu, 'String');
    for i = numel(handles.ROILegend):-1:1
        try
            handles.ROILegend(i) = [];
        end
    end
    for i = 1:numel(handles.ROIs)
        if i > 1
            box     = get( handles.ROILegend{i-1},'Extent' );
%             yPos    = box(2)-box(4);
            yPos    = box(2)-handles.ROILegendHeight;
        else
            yPos    = handles.yMax;
        end
        handles.ROILegend{i} = text( 'Position', [1,yPos,0], 'String', legendStr{i}, 'FontSize', 12,'Color',...
            handles.ROIColor{mod(i,numel(handles.ROIColor))+1},...
            'HorizontalAlignment','Left','VerticalAlignment','bottom', ...
            'interpreter','none');
    end
    
    % Left-right info
    text( 'String','R', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.005, 0.5, 0],'FontSize',10);
    text( 'String','L', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.980, 0.5, 0],'FontSize',10);
    text( 'String','A', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.985, 0],'FontSize',10);
    text( 'String','P', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.02, 0],'FontSize',10);
    
    % Draw data
    handles       = drawImagePlusROI(handles);
    guidata(hObject, handles);
    if nargout > 0 
        varargout = {handles};
    end
end

% --- Executes on button press in sagPushbutton.
function sagPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to sagPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    handles.proj = 'sag';
    set(handles.projText,'String','Sagittal (XZ slice)'); 
    %Delete line profile if any
    delete( findall( handles.axes1,'type','line') );
    set( handles.textLineProfWidth,'Enable', 'off' );
    set( handles.sliderLineProfWidth,'Enable', 'off' );
    set( handles.textLineProfText,'Enable', 'off' );
    set( handles.textLineProfWidth,'String', [num2str( get(handles.sliderLineProfWidth,'Value') )] );
    
    handles.sliderValueSpace     = round( get(handles.sliceSlider,'Value') );
    if handles.sliderValueSpace > size(handles.data,2)
        handles.sliderValueSpace = size(handles.data,2);
    end
    
    % Set slider values
    % If the data is 3D...
    if size(handles.data,2)>1
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,2)-1)) 10/((size(handles.data,2)-1))]);
    %If it isn't...
    else
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,2))) 10/((size(handles.data,2)))]);
    end
    set(handles.sliceSlider,'Max',size(handles.data,2));
    set(handles.textTotalSlice,'String',[ '/ ', num2str( size(handles.data,2) )] );
    set(handles.sliceSlider,'value', handles.sliderValueSpace);
    sliceSlider_Callback( handles.sliceSlider,[],handles);
    
    % Inactivate sliders if max = min
    if size(handles.data,2)==1
        handles.sliceSlider.Enable = 'off';
    else
        handles.sliceSlider.Enable = 'on';
    end

    % Update text info at top
    updateMinMaxText(handles);
    
    % What is shown on the image display
    handles.xMax = size(handles.data,1);
    handles.yMax = size(handles.data,3);

    % Choose axes to draw in
    axes(handles.axes1);

    % Plot image
    imshow( rot90( squeeze(handles.data(:,round(handles.sliderValueSpace),:,round(handles.sliderValueTime))) ),[handles.lowLevel handles.highLevel]);
    daspect( [ handles.voxelSize(3) handles.voxelSize(1) 1] );
    colormaps_Callback(handles.colormaps, [], handles)
    axis off

    % Text in axes
    legendStr = get( handles.ROINumberMenu, 'String');
    for i = numel(handles.ROILegend):-1:1
        try
            handles.ROILegend(i) = [];
        end
    end
    for i = 1:numel(handles.ROIs)
        if i > 1
            box     = get( handles.ROILegend{i-1},'Extent' );
            yPos    = box(2)-handles.ROILegendHeight;
        else
            yPos    = handles.yMax;
        end
        handles.ROILegend{i} = text( 'Position', [1,yPos,0], 'String', legendStr{i}, 'FontSize', 12,'Color',...
            handles.ROIColor{mod(i,numel(handles.ROIColor))+1},...
            'HorizontalAlignment','Left','VerticalAlignment','bottom', ...
            'interpreter','none');
    end
    
    % Left-right info
    text( 'String','A', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.005, 0.5, 0],'FontSize',10);
    text( 'String','P', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.980, 0.5, 0],'FontSize',10);
    text( 'String','S', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.985, 0],'FontSize',10);
    text( 'String','I', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.02, 0],'FontSize',10);
    
    % Draw data
    handles = drawImagePlusROI(handles);
    guidata(hObject, handles);
end

% --- Executes on button press in corPushbutton.
function corPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to corPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    handles.proj = 'cor'; 
    set(handles.projText,'String','Coronal (YZ slice)'); 

    %Delete line profile if any
    delete( findall( handles.axes1,'type','line') );
    set( handles.textLineProfWidth,'Enable', 'off' );
    set( handles.sliderLineProfWidth,'Enable', 'off' );
    set( handles.textLineProfText,'Enable', 'off' );
    set( handles.textLineProfWidth,'String', [num2str( get(handles.sliderLineProfWidth,'Value') )] );
    
    % Set slider max value
    handles.sliderValueSpace     = round( get(handles.sliceSlider,'Value') );
    if handles.sliderValueSpace > size(handles.data,1)
        handles.sliderValueSpace = size(handles.data,1);
    end
    % Set slider values
    % If the data is 3D...
    if size(handles.data,1)>1
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,1)-1)) 10/((size(handles.data,1)-1))]);
    %If it isn't...
    else
        set(handles.sliceSlider,'SliderStep',[1/((size(handles.data,1))) 10/((size(handles.data,1)))]);
    end
    set(handles.sliceSlider,'Max',size(handles.data,1));
    set(handles.textTotalSlice,'String',[ '/ ', num2str( size(handles.data,1) )] );
    set(handles.sliceSlider,'value', handles.sliderValueSpace);
    sliceSlider_Callback( handles.sliceSlider,[],handles);
    
    % Inactivate sliders if max = min
    if size(handles.data,1)==1
        handles.sliceSlider.Enable = 'off';
    else
        handles.sliceSlider.Enable = 'on';
    end

    % Update text info at top
    updateMinMaxText(handles);
    
    handles.xMax = size(handles.data,2);
    handles.yMax = size(handles.data,3);

    % Choose axes to draw in
    axes(handles.axes1);
    
    % Plot image
    imshow( rot90( squeeze(handles.data(round(handles.sliderValueSpace),:,:,round(handles.sliderValueTime))) ),[handles.lowLevel handles.highLevel]);
    daspect( [ handles.voxelSize(3) handles.voxelSize(2) 1] );
    colormaps_Callback(handles.colormaps, [], handles)
    axis off

    % Text in axes
    legendStr = get( handles.ROINumberMenu, 'String');
    for i = numel(handles.ROILegend):-1:1
        try
            handles.ROILegend(i) = [];
        end
    end
    for i = 1:numel(handles.ROIs)
        if i > 1
            box     = get( handles.ROILegend{i-1},'Extent' );
            yPos    = box(2)-handles.ROILegendHeight;
        else
            yPos    = handles.yMax;
        end
        handles.ROILegend{i} = text( 'Position', [1,yPos,0], 'String', legendStr{i}, 'FontSize', 12,'Color',...
            handles.ROIColor{mod(i,numel(handles.ROIColor))+1},...
            'HorizontalAlignment','Left','VerticalAlignment','bottom',...
            'interpreter','none');
    end
    
    % Left-right info
    text( 'String','R', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.005, 0.5, 0],'FontSize',10);
    text( 'String','L', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.980, 0.5, 0],'FontSize',10);
    text( 'String','S', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.985, 0],'FontSize',10);
    text( 'String','I', 'Parent',handles.axes1,'Color',[0.8 0.8 0.8],...
        'Units','Normalized','Position',[0.5, 0.02, 0],'FontSize',10);
    
    % Draw data
    handles = drawImagePlusROI(handles);
    guidata(hObject, handles);
end

%==========================================================================
%                               Sliders
%==========================================================================

%============
% Sliceslider
%============

% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
    % hObject    handle to sliceSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Obtains the slider value from the slider component
    handles.sliderValueSpace = round(get(handles.sliceSlider,'Value'));
    %Calculate corresponding slice in length units
    switch lower(handles.proj)
        case 'tra'
            sliceUnit = [ num2str( handles.slicePosition(handles.sliderValueSpace) ) ' ' handles.dimUnits{1} ];
        case 'sag'
            sliceUnit = [ num2str( (handles.sliderValueSpace - size(handles.data,2)/2 - ...
                1/2 ) * handles.voxelSize(2) ) ' ' handles.dimUnits{1} ];
        case 'cor'
            sliceUnit = [ num2str( (handles.sliderValueSpace - size(handles.data,1)/2 - ...
                1/2 ) * handles.voxelSize(1) ) ' ' handles.dimUnits{1} ];
    end  
    set(handles.sliceUnitText,'String', sliceUnit );
    
    %puts the slider value into the edit text component
    set(handles.sliceNumber,'String', num2str(handles.sliderValueSpace));

    % Update text info at top
    updateMinMaxText(handles);
    
    % Update line profiles (if any)
    if ~isempty( findall(handles.axes1,'type','line') )
        colCoord = get( handles.lineProfileCenter, 'XData');
        rowCoord = get( handles.lineProfileCenter, 'YData');
        width    = get( handles.sliderLineProfWidth, 'Value');
        handles  = drawLineProfileInAxis(handles,colCoord(1),rowCoord(1),width);
    end
    
    handles = drawImagePlusROI(handles);

    guidata(hObject, handles);
end


%============
% Timeslider
%============

% --- Executes on slider movement.
function timeSlider_Callback(hObject, eventdata, handles)
    % hObject    handle to timeSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Obtains the slider value from the slider component
    handles.sliderValueTime = round(get(handles.timeSlider,'Value'));
    
    %Calculate corresponding tim in time units
    timeUnit = [ num2str( handles.timeMidFrame( handles.sliderValueTime ) ) ' ' handles.dimUnits{2} ];
    set(handles.timeUnitText,'String', timeUnit );
    
    %puts the slider value into the edit text component
    set(handles.timeNumber,'String', num2str(handles.sliderValueTime));
    
    % Update text info at top
    updateMinMaxText(handles);
    
    % Update line profiles (if any)
    if ~isempty( findall(handles.axes1,'type','line') )
        colCoord = get( handles.lineProfileCenter, 'XData');
        rowCoord = get( handles.lineProfileCenter, 'YData');
        width    = get( handles.sliderLineProfWidth, 'Value');
        handles  = drawLineProfileInAxis(handles,colCoord(1),rowCoord(1),width);
    end
    
    handles = drawImagePlusROI(handles);
    
    %If ROI is chosen, update ROI statistics
    if strcmp(get(handles.ROINumberMenu,'Enable'),'on')
        handles = calculateROIStatistics(handles);
    end
    
    
    guidata(hObject, handles);
end


%============
% lowslider
%============

% --- Executes on slider movement.
function lowSlider_Callback(hObject, eventdata, handles)
    % hObject    handle to lowSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Obtains the slider value from the slider component
    sliderValueLow = (get(handles.lowSlider,'Value'));
    sliderValueHigh = (get(handles.highSlider,'Value'));

    % Checks that the lower threshold is <= upper threshold
    if sliderValueLow >= sliderValueHigh
        step = get(handles.highSlider,'sliderStep');
        handles.lowLevel = sliderValueHigh - step(1);
        set(handles.lowSlider,'Value',handles.lowLevel);
    else
        handles.lowLevel = sliderValueLow;   
    end

    %puts the new slider value into the edit text component
    set(handles.lowNumber,'String',num2str(handles.lowLevel))

    %Update Graph
    set( handles.axes1,'CLim',[handles.lowLevel handles.highLevel]);
    handles = drawImagePlusROI(handles);
    guidata(hObject, handles);
end


%============
% highslider
%============
% --- Executes on slider movement.
function highSlider_Callback(hObject, eventdata, handles)
    % hObject    handle to highSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Obtains the slider value from the slider component
    sliderValueLow = (get(handles.lowSlider,'Value'));
    sliderValueHigh = (get(handles.highSlider,'Value'));

    % Checks that the lower threshold is <= upper threshold
    if sliderValueLow >= sliderValueHigh
        step = get(handles.highSlider,'sliderStep');
        handles.highLevel = sliderValueLow + step(1);
        set(handles.highSlider,'Value',handles.highLevel);
    else
        handles.highLevel = sliderValueHigh;   
    end

    %puts the slider value into the edit text component
    set(handles.highNumber,'String',num2str(handles.highLevel))
    
    %Update display
    set( handles.axes1,'CLim',[handles.lowLevel handles.highLevel]);
    handles = drawImagePlusROI(handles);
    
    guidata(hObject, handles);
end


%==========================================================================
%                                   Text
%==========================================================================

%============
% Slicenumber
%============
function sliceNumber_Callback(hObject, eventdata, handles)
    % hObject    handle to sliceNumber (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of sliceNumber as text
    %        str2double(get(hObject,'String')) returns contents of sliceNumber as a double

    %get the string for the editText component
    handles.sliderValueSpace = str2num(get(handles.sliceNumber,'String'));

    %if user inputs something is not a number, or if the input is less than the
    %min or greater than the max value these statements catch it.
    if (isempty(handles.sliderValueSpace) || handles.sliderValueSpace < get(handles.sliceSlider,'Min'))
        set(handles.sliceSlider,'Value',get(handles.sliceSlider,'Min'));
    elseif handles.sliderValueSpace > get(handles.sliceSlider,'Max')
        set((handles.sliceSlider),'Value',get(handles.sliceSlider,'Max'));
    else
        set(handles.sliceSlider,'Value',handles.sliderValueSpace);
    end
    
    % Update text info at top
    updateMinMaxText(handles);
    
    guidata(hObject, handles);
    sliceSlider_Callback( handles.sliceSlider,[],handles);
end


%============
% Timenumber
%============
function timeNumber_Callback(hObject, eventdata, handles)
    % hObject    handle to timeNumber (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    handles.sliderValueTime = str2num(get(handles.timeNumber,'String'));

    %if user inputs something is not a number, or if the input is less than 0
    %or greater than 100, then the slider value defaults to 0
    if (isempty(handles.sliderValueTime) || handles.sliderValueTime < get(handles.timeSlider,'Min'))
        set(handles.timeSlider,'Value',get(handles.timeSlider,'Min'));
    elseif handles.sliderValueTime > get(handles.timeSlider,'Max')
        set((handles.timeSlider),'Value',get(handles.timeSlider,'Max'));
    else
        set(handles.timeSlider,'Value',handles.sliderValueTime);
    end
    
    % Update text info at top
    updateMinMaxText(handles);
    
    guidata(hObject,handles)
    timeSlider_Callback( handles.timeSlider,[],handles);
end


%============
% Lowernumber
%============

function lowNumber_Callback(hObject, eventdata, handles)
    % hObject    handle to lowNumber (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    sliderValueLow = (str2num(get(handles.lowNumber,'String')));
    sliderValueHigh = (get(handles.highSlider,'Value'));

    %if user inputs something is not a number, or if the input is less than 0
    %or greater than 100, then the slider value defaults to 0
    if (isempty(sliderValueLow) || sliderValueLow < get(handles.lowSlider,'Min'))
        handles.lowLevel = get(handles.lowSlider,'Min');
    elseif sliderValueLow >= sliderValueHigh
        step = get(handles.lowSlider,'sliderStep');
        handles.lowLevel = sliderValueHigh - step(1);   
    else
        handles.lowLevel = sliderValueLow;   
    end
    set(handles.lowSlider,'Value',handles.lowLevel);
    set(handles.lowNumber,'String',num2str(handles.lowLevel))

    %Update display
    set( handles.axes1,'CLim',[handles.lowLevel handles.highLevel]);
    handles = drawImagePlusROI(handles);

    guidata(hObject,handles)
end


%============
% Uppernumber
%============

function highNumber_Callback(hObject, eventdata, handles)
    % hObject    handle to highNumber (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    sliderValueLow  = (get(handles.lowSlider,'Value'));
    sliderValueHigh = (str2num(get(handles.highNumber,'String'))); 

    %if user inputs something is not a number, or if the input is less than 0
    %or greater than max of image, then the slider value defaults to max
    if (isempty(sliderValueHigh) || sliderValueHigh > get(handles.highSlider,'Max'))
        handles.highLevel = get(handles.highSlider,'Max');
    elseif (sliderValueLow >= sliderValueHigh)
        step = get(handles.highSlider,'sliderStep');
        handles.highLevel = sliderValueLow + step(1);
    else
        handles.highLevel = sliderValueHigh;
    end
    set(handles.highSlider,'Value',handles.highLevel);
    set(handles.highNumber,'String',num2str(handles.highLevel))
    
    %Update display
    set( handles.axes1,'CLim',[handles.lowLevel handles.highLevel]);
    handles = drawImagePlusROI(handles);

    guidata(hObject,handles)
end


%============
% Colormap
%============
function colormaps_Callback(hObject, eventdata, handles)
    % hObject    handle to colormaps (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    names = get(hObject,'String');
    colormap( handles.figure1, names{get(hObject,'Value')} );
    if isfield(handles,'colorbar')
        tmp = get(handles.colormaps,'String');
        colormap( handles.colorbar, tmp{get(handles.colormaps,'Value')} );
    end
    mycmap = get(handles.figure1,'Colormap');
    if (get(handles.invertCheckbox,'Value') == 1)
%         set(handles.figure1,'Colormap',flipud(mycmap))
%         set(handles.colorbar,'Colormap',flipud(mycmap))
        colormap( handles.figure1, flipud(mycmap) );
        colormap( handles.colorbar, flipud(mycmap) );
    end
    

    guidata(hObject,handles)
end


% --- Executes on button press in TACTpushbutton.
function TACTpushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to TACTpushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    %If no ROI chosen, do nothing
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end
   
    % Screen output
    disp('Generating TACs...')
    
    % Calculate ROI values
    for i = 1:numel(handles.ROIs)
        
        % ROIed image
        ROIindex        = find( handles.ROIs{i} > 0 );       
        for j = 1:size( handles.data,4 )
            image3D     = handles.data(:,:,:,j);
            TAC(j,i)    = mean( image3D(ROIindex) );
        end
        
        % Time vector
        time(:,i)   = handles.timeMidFrame;
        
        % Check if you should scale to framelength
        if ( get(handles.TACTScale2TimeCheckbox,'Value') == 1 )
            TAC(:,i) = TAC(:,i) ./ handles.frameLength;
            yLabel   = ['Intensity  (' handles.imageUnit '/' handles.dimUnits{2} ')'];
        else
            yLabel   = ['Intensity  (' handles.imageUnit ')'];
        end
        
        % Plot figure in existing or new window
        fh = findobj('type','figure','name','TAC');
        if isempty( fh )
            fh = figure('name','TAC','numbertitle','off'); 
            plot( time(:,i), TAC(:,i), '-s','Color','k','LineWidth',1,'MarkerSize',5,'MarkerFaceColor',...
                handles.ROIColor{mod(i,numel(handles.ROIColor))+1},'MarkerEdgeColor',[0 0 0]);
            hold on;
        else
            figure(fh);
            plot( time(:,i), TAC(:,i), '-s','Color','k','LineWidth',1,'MarkerSize',5,'MarkerFaceColor',...
                handles.ROIColor{mod(i,numel(handles.ROIColor))+1},'MarkerEdgeColor',[0 0 0]);
            hold on;
        end

    end
    figure(fh); hold off;
    xlabel(['Time  (' handles.dimUnits{2} ')'],'FontSize',12), ylabel(yLabel,'FontSize',12), title('TAC','FontSize',12);
    legendStr = get( handles.ROINumberMenu,'String');
    legend( legendStr,'Interpreter','none' );
    
    % Set focus back to ROIViewer4D
    set(handles.figure1,'Visible','on')
    
    % Check if you should export to workspace
    if ( get(handles.TACExportCheckbox,'Value') == 1 )
        for i = 1:size(TAC,2)
            if i == 1
                assignin('base',['TAC_t'],time(:,i));
            end
            assignin('base',['TAC_' legendStr{i}],TAC(:,i));
        end
    end
    
    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles)
end


% --- Executes on button press in hideROIcheckbox.
function hideROIcheckbox_Callback(hObject, eventdata, handles)
    % hObject    handle to hideROIcheckbox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    handles = drawImagePlusROI(handles);
    guidata(hObject,handles)
end


% --------------------------------------------------------------------
% ROI drawing functions
% --------------------------------------------------------------------
function clickImageFunction(src, eventdata)
    handles = guidata(src);

    %If right click, do nothing
    if strcmp(get(src,'SelectionType'),'alt')
        return
    end

    % If clicking outside image, do nothing
    coordinates = round( get(handles.axes1,'CurrentPoint') );
    if ( coordinates(1,1) < handles.xMin ) || ( coordinates(1,1) > handles.xMax ) ...
            || ( coordinates(1,2) < handles.yMin ) || ( coordinates(1,2) > handles.yMax )
        return
    end
    
    % Line profile coordinates
    if ~isempty( findall( handles.axes1,'type','line') )
        lineProfX = get( handles.lineProfileCenter,'XData');
        lineProfY = get( handles.lineProfileCenter,'YData');
    else
        lineProfX = [-1 -1];
        lineProfY = [-1 -1];
    end
    
    % If clicking on a line profile...
    if ( coordinates(1,1) >= lineProfX(1) ) && ( coordinates(1,1) <= lineProfX(2) ) ...
            && ( coordinates(1,2) >= lineProfY(1) ) && ( coordinates(1,2) <= lineProfY(2) )
        handles = clickLineProfile(src, eventdata);
        guidata(src,handles);
        return
    elseif strcmp(get(handles.ROINumberMenu,'Enable'),'off') %If no ROI selected, do nothing
        return
    end
    
    %If clicking to draw ROI
    set(src,'windowButtonUpFcn',@unClickImageFunction)
    set(src,'windowButtonMotionFcn',{@ROIMovePoint})

    handles = updateROI(handles);
    handles = drawImagePlusROI(handles);

    guidata(src,handles);
end


function handles = clickLineProfile(src, eventdata)
    handles = guidata(src);
    
    set(src,'windowButtonMotionFcn',{@lineProfileMovePoint})
    set(src,'windowButtonUpFcn',@unClickLineProfile)

    handles = updateLineProfile(handles);
    
    guidata(src,handles);    
end


function ROIMovePoint(src, eventdata)
    handles = guidata(src);

    set(src,'windowButtonDownFcn',@clickImageFunction)
    set(src,'windowButtonUpFcn',@unClickImageFunction)

    coordinates = round( get(handles.axes1,'CurrentPoint') );
    if ( coordinates(1,1) >= handles.xMin ) && ( coordinates(1,1) <= handles.xMax ) ...
            && ( coordinates(1,2) >= handles.yMin ) && ( coordinates(1,2) <= handles.yMax )
        
        t       = str2num( get(handles.timeNumber,'String') );
        slice   = str2num( get(handles.sliceNumber,'String') );
        x       = coordinates(1,1); % horizontal in displayed image
        y       = coordinates(1,2); % vertical in displayed image
        if handles.proj == 'tra'
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(1) size( handles.data,2) size( handles.data,1) ];
            value     = handles.data(y,x,slice,t);
        elseif handles.proj == 'sag'
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(1) handles.voxelSize(3) size( handles.data,1) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(x,slice,y,t);
        elseif handles.proj == 'cor'
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(3) size( handles.data,2) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(slice,x,y,t);
        end
        
        xUnit       = ( x - pixelSize(3)/2 - 1/2 ) * pixelSize(1);
        yUnit       = ( y - pixelSize(4)/2 - 1/2 ) * pixelSize(2);
        if ( abs(value) >= 1e4 ) %|| ( value <= 1e-3 )
                string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value,'%.3e') ];
        else
                string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value) ];
        end
        stringUnit  = ['(' num2str(xUnit,'%.2f') ', ' num2str(yUnit,'%.2f') ') ' handles.dimUnits{1} ];
        set( handles.textPixelInfo,'String',string );
        set( handles.textPixelInfoUnit,'String',stringUnit );
    end
    
    % Update ROI
    handles = updateROI(handles);

    guidata(src,handles);
end


function lineProfileMovePoint(src, eventdata)
    handles = guidata(src);
    
    set(src,'windowButtonDownFcn',@clickImageFunction)
    set(src,'windowButtonUpFcn',@unClickLineProfile)

    coordinates = round( get(handles.axes1,'CurrentPoint') );
    if ( coordinates(1,1) >= handles.xMin ) && ( coordinates(1,1) <= handles.xMax ) ...
            && ( coordinates(1,2) >= handles.yMin ) && ( coordinates(1,2) <= handles.yMax )
        
        t       = str2num( get(handles.timeNumber,'String') );
        slice   = str2num( get(handles.sliceNumber,'String') );
        x       = coordinates(1,1); % horizontal in displayed image
        y       = coordinates(1,2); % vertical in displayed image
        if handles.proj == 'tra'
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(1) size( handles.data,2) size( handles.data,1) ];
            value     = handles.data(y,x,slice,t);
        elseif handles.proj == 'sag'
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(1) handles.voxelSize(3) size( handles.data,1) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(x,slice,y,t);
        elseif handles.proj == 'cor'
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(3) size( handles.data,2) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(slice,x,y,t);

        end
        
        xUnit       = ( x - pixelSize(3)/2 - 1/2 ) * pixelSize(1);
        yUnit       = ( y - pixelSize(4)/2 - 1/2 ) * pixelSize(2);
        if ( abs(value) >= 1e4 ) %|| ( value <= 1e-3 )
                string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value,'%.3e') ];
        else
                string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value) ];
        end
        stringUnit  = ['(' num2str(xUnit,'%.2f') ', ' num2str(yUnit,'%.2f') ') ' handles.dimUnits{1} ];
        set( handles.textPixelInfo,'String',string );
        set( handles.textPixelInfoUnit,'String',stringUnit );
    end
    
    % Update line profile
    handles = updateLineProfile(handles);

    guidata(src,handles);
end


function handles = updateROI(handles)

    %Current pointer location
    coordinate = get(handles.axes1,'CurrentPoint');
    x          = round( coordinate(1,1) ); %horizontal
    y          = round( coordinate(1,2) ); %vertical
    
    %Current settings
    ROINumber = get(handles.ROINumberMenu,'Value');
    sliceNo   = round( get(handles.sliceSlider,'Value') );
    imageSize = size( get(findall(handles.axes1,'type','image'),'CData') );
    
    % Brush size around pointer location
    r         = str2num( get(handles.BrushSize,'String') );
        
    % If set to "Brush" or "Eraser"
    if ( get(handles.ROIBrushRadiobutton,'Value') == 1 )
        roiValue = 1; % Brush
    else
        roiValue = -1; % Eraser
    end
                
    
    % Check type of brush (circle or square)
    if get(handles.radiobuttonCircleROI,'Value') == 1
        % Circular brush
        [Y,X] = meshgrid( -(x-1):(imageSize(2)-x), -(y-1):(imageSize(1)-y) );
        mask  = ((X.^2+Y.^2) <= r^2);
    else
        % Square brush
        mask  = zeros( imageSize(1),imageSize(2) );
        xMinInd = round(x-(r-1)/2);
        xMaxInd = round(x+(r-1)/2);
        yMinInd = round(y-(r-1)/2);
        yMaxInd = round(y+(r-1)/2);
        if xMaxInd > handles.xMax
            xMaxInd = handles.xMax;
        end
        if xMinInd < handles.xMin
            xMinInd = handles.xMin;
        end
        if yMaxInd > handles.yMax
            yMaxInd = handles.yMax;
        end
        if yMinInd < handles.yMin
            yMinInd = handles.yMin;
        end
        mask(yMinInd:yMaxInd,xMinInd:xMaxInd)  = 1;
    end
            
    %Data to put in ROI
    switch lower(handles.proj)
        case 'tra'
            %ROI
            tempData                 = squeeze( handles.ROIs{ROINumber}(:,:,sliceNo) );
            tempData                 = roiValue*mask + tempData;
            tempData( tempData > 1 ) = 1;
            tempData( tempData < 0 ) = 0;
            handles.ROIs{ROINumber}(:,:,sliceNo) = tempData;
        case 'cor'
            %ROI, displayed image rotated 90 degrees compared to ROI array
            tempData                 = squeeze( handles.ROIs{ROINumber}(sliceNo,:,:) );
            tempData                 = roiValue*rot90(mask,3) + tempData;
            tempData( tempData > 1 ) = 1;
            tempData( tempData < 0 ) = 0;
            handles.ROIs{ROINumber}(sliceNo,:,:) = tempData;
        case 'sag'
            %Displayed image rotated 90 degrees compared to ROI array
            tempData                 = squeeze( handles.ROIs{ROINumber}(:,sliceNo,:) );
            tempData                 = roiValue*rot90(mask,3) + tempData;
            tempData( tempData > 1 ) = 1;
            tempData( tempData < 0 ) = 0;
%             handles.ROIs{ROINumber}(:,sliceNo,:,timeNo) = rot90( tempData );
            handles.ROIs{ROINumber}(:,sliceNo,:) = tempData;
    end  
       
    handles = drawImagePlusROI(handles);
end


function handles = updateLineProfile(handles)

    %Current pointer location
    coordinate      = round( get(handles.axes1,'CurrentPoint') );
    colCoord        = coordinate(1,1);
    rowCoord        = coordinate(1,2);
    width           = round( get( handles.sliderLineProfWidth,'Value' ) );
 
    % Can't go outside image borders
    if mod(width,2) == 0
        if ( colCoord - width/2 < handles.xMin )
            colCoord = handles.xMin + width/2 - 1;
        elseif ( colCoord + width/2 > handles.xMax )
            colCoord = handles.xMax - width/2;
        end
        if ( rowCoord - width/2 < handles.yMin )
            rowCoord = handles.yMin + width/2 - 1;
        elseif ( rowCoord + width/2 > handles.yMax )
            rowCoord = handles.yMax - width/2;
        end
    else
        if ( colCoord - width/2 < handles.xMin )
            colCoord = handles.xMin + floor(width/2);
        elseif ( colCoord + width/2 > handles.xMax )
            colCoord = handles.xMax - floor(width/2);
        end
        if ( rowCoord - width/2 < handles.yMin )
            rowCoord = handles.yMin + floor(width/2);
        elseif ( rowCoord + width/2 > handles.yMax )
            rowCoord = handles.yMax - floor(width/2);
        end
    end
    handles = drawLineProfileInAxis(handles,colCoord,rowCoord,width);
    
end

function handles = drawLineProfileInAxis(handles,colCoord,rowCoord,width)
       
    switch handles.lineProfileDirection
        case 'Vertical'
            if mod(width,2) == 0
                set( handles.lineProfile,'Position',[colCoord-width/2+1 0.5 width-1 handles.yMax-0.5]);
            else
                set( handles.lineProfile,'Position',[colCoord-width/2 0.5 width handles.yMax-0.5]);
            end
            set( handles.lineProfileCenter,'XData',[colCoord colCoord]);
            set( handles.lineProfileCenter,'YData',[0 handles.yMax+1]);
        case 'Horizontal'
            if mod(width,2) == 0
                set( handles.lineProfile,'Position',[0.5 rowCoord-width/2+1 handles.xMax-0.5 width-1]);
            else
                set( handles.lineProfile,'Position',[0.5 rowCoord-width/2 handles.xMax-0.5 width]);
            end
            set( handles.lineProfileCenter,'XData',[0 handles.xMax+1]);
            set( handles.lineProfileCenter,'YData',[rowCoord rowCoord]);
    end
      
    % Plot line profile
    sliceNo = str2num( get(handles.sliceNumber,'String') );
    timeNo  = str2num( get(handles.timeNumber,'String') );
    
    switch lower(handles.proj)
        case 'tra'
            image = squeeze( handles.data(:,:,sliceNo,timeNo) );
        case 'cor'
            image = rot90( squeeze( handles.data(sliceNo,:,:,timeNo) ),1 );
        case 'sag'
            image = rot90( squeeze( handles.data(:,sliceNo,:,timeNo) ),1 );
    end
    
    switch handles.lineProfileDirection
        case 'Vertical'
            % Check if even or uneven line profile width
            if mod(width,2) == 0
                index = colCoord-width/2+1:colCoord+width/2;
            else
                index = colCoord-floor(width/2):colCoord+floor(width/2);
            end
        case 'Horizontal'
            % Check if even or uneven line profile width
            if mod(width,2) == 0
                index = rowCoord-width/2+1:rowCoord+width/2;
            else
                index = rowCoord-floor(width/2):rowCoord+floor(width/2);
            end
    end
    
    handles = plotLineProfile(handles,image,index);
    
end


function unClickImageFunction(src, eventdata)
    handles = guidata(src);

    % It right click, nothing happens
    if strcmp(get(src,'SelectionType'),'alt')
        return
    end
    %If no ROI selected, do nothing
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    set(src,'windowButtonMotionFcn',@(hObject,eventdata)ROIViewer4D('figure1_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject)))
    set(src,'KeyPressFcn','')

    handles = calculateROIStatistics(handles);
    guidata(src, handles);
end


function unClickLineProfile(src, eventdata)
    handles = guidata(src);
    
    set(src,'windowButtonUpFcn',@unClickImageFunction)
    set(handles.figure1,'WindowButtonDownFcn',@clickImageFunction); 
    set(src,'windowButtonMotionFcn',@(hObject,eventdata)ROIViewer4D('figure1_WindowButtonMotionFcn',hObject,eventdata,guidata(hObject)))
    set(src,'KeyPressFcn','');
           
    guidata(src, handles);
end


function handles = calculateROIStatistics(handles)

    % Calculate ROI statistics
    ROINumber   = get(handles.ROINumberMenu,'Value');
    timeNo      = str2num( get( handles.timeNumber,'String') );
    ROIindex    = find( handles.ROIs{ROINumber} > 0 );
    try 
        if ~isempty( ROIindex )
            ROIedTemp   = handles.data(:,:,:,timeNo);
            ROIedImage  = ROIedTemp(ROIindex);
        else
            ROIedImage  = 0;
        end        
    catch
        ROIedImage  = handles.data;
    end
    
    ROIVolume       = numel(ROIindex) * handles.voxelSize(1) * handles.voxelSize(2) * handles.voxelSize(3);
    if strcmp( handles.dimUnits(1), 'mm') == 1
        ROIVolume   = ROIVolume/1000; % To get it in [cm^3]
        voxelString = [ num2str(numel(ROIindex)) '   (' num2str(ROIVolume,'%.2f') ' cm^3)' ];
    else
        voxelString = [ num2str(numel(ROIindex)) '   (' num2str(ROIVolume,'%.2f') ' ' handles.dimUnits{1} '^3)' ];
    end

    if mean( ROIedImage ) < 1e4
        set( handles.textROIMean,'String',num2str(mean( ROIedImage )) );
    else
        set( handles.textROIMean,'String',num2str(mean( ROIedImage ),'%.3e') );
    end
    if std( ROIedImage ) < 1e4
        if mean(ROIedImage) == 0
            set( handles.textROIStd,'String',[num2str(std( ROIedImage )) ...
            '   (' num2str(100*std(ROIedImage),'%.2f') '%)'] );
        else
            set( handles.textROIStd,'String',[num2str(std( ROIedImage )) ...
            '   (' num2str(100*std(ROIedImage)/mean(ROIedImage),'%.2f') '%)'] );
        end
    else
        if mean(ROIedImage) == 0
            set( handles.textROIStd,'String',[num2str(std( ROIedImage ),'%.3e') ...
            '   (' num2str(100*std(ROIedImage),'%.2f') '%)'] );
        else
            set( handles.textROIStd,'String',[num2str(std( ROIedImage ),'%.3e') ...
            '   (' num2str(100*std(ROIedImage)/mean(ROIedImage),'%.2f') '%)'] );
        end
    end
    if min( ROIedImage ) < 1e4
        set( handles.textROIMin,'String',num2str(min( ROIedImage )) );
    else
        set( handles.textROIMin,'String',num2str(min( ROIedImage ),'%.3e') );
    end
    if max( ROIedImage ) < 1e4
        set( handles.textROIMax,'String',num2str(max( ROIedImage )) );
    else
        set( handles.textROIMax,'String',num2str(max( ROIedImage ),'%.3e') );
    end
    if sum( ROIedImage(:) ) < 1e4
        set( handles.textROISum,'String',num2str(sum( ROIedImage(:) )) );
    else
        set( handles.textROISum,'String',num2str(sum( ROIedImage(:) ),'%.3e') );
    end
    set( handles.textROINoVoxels,'String', voxelString );

end


function figScroll(src, eventdata)
    % guidata(hObject, handles);
end


function handles = drawImagePlusROI(handles)
    % Function that draws the image and chosen ROI

    % Current values
    ROINumber   = get(handles.ROINumberMenu,'Value');
    sliceNo     = round( get(handles.sliceSlider,'Value') );
    timeNo      = round( get(handles.timeSlider,'Value') );

    % Extract 2D image (currently shown) and corresponding 2D ROI
    if strcmp(handles.proj,'tra')
        current2DImage   = squeeze( handles.data(:,:,sliceNo,timeNo) );
        current2DROIFlat = zeros(size(current2DImage));
        try
            for i = 1:numel( handles.ROIs )
                current2DROI{i}  = squeeze( handles.ROIs{i}(:,:,sliceNo) );
                current2DROIFlat = current2DROIFlat + current2DROI{i};
            end
        catch
            current2DROI = zeros( size(current2DImage) );
        end
        
    elseif strcmp(handles.proj,'cor')
        % Rotate displayed image 90 degreees
        current2DImage   = rot90( squeeze( handles.data(sliceNo,:,:,timeNo) ) );
        current2DROIFlat = zeros(size(current2DImage));
        try
            for i = 1:numel( handles.ROIs )
                current2DROI{i}  = rot90( squeeze( handles.ROIs{i}(sliceNo,:,:) ) );
                current2DROIFlat = current2DROIFlat + current2DROI{i};
            end
        catch
            current2DROI = zeros( size(current2DImage) );
        end
        
    elseif strcmp(handles.proj,'sag')
        % Rotate displayed image 90 degreees
        current2DImage   = rot90( squeeze( handles.data(:,sliceNo,:,timeNo) ) );
        current2DROIFlat = zeros(size(current2DImage));
        try
            for i = 1:numel( handles.ROIs )
                current2DROI{i}  = rot90( squeeze( handles.ROIs{i}(:,sliceNo,:) ) );
                current2DROIFlat = current2DROIFlat + current2DROI{i};
            end
        catch
            current2DROI = zeros( size(current2DImage) );
        end
        
    end

    % If "Hide ROI" is checked, use ROI of only zeros
    if ( get(handles.hideROIcheckbox,'Value') == 1 )
        current2DROIFlat = zeros( size(current2DImage) );
        set( findall( handles.axes1,'Tag','Contour'),'Visible','off' );
    end

    % If "Show only ROI" is checked, use image of only zeros
    if ( get(handles.showOnlyROICheckbox,'Value') == 1 )
        current2DImage = zeros( size(current2DImage) );
        %Hide line profile
        if ~isempty( findall( handles.axes1,'type','line') );
            set( findall( handles.axes1,'type','line'), 'Visible', 'off')
        end
    else
        %Show line profile
        if ~isempty( findall( handles.axes1,'type','line') );
            set( findall( handles.axes1,'type','line'), 'Visible', 'on')
        end
    end

    % Image to be displayed; image + ROI
    if max(current2DImage(:)) ~= 0
        combinedImage = current2DImage + (handles.ROIwindow*max(current2DImage(:))) * current2DROIFlat;
    else
        combinedImage = current2DImage + current2DROIFlat;
    end
    % If scale image to slice 
    if ( get( handles.Scale2SliceRadioButton,'Value') == 1 )
        % Image max and min value
        if handles.proj == 'tra'
            sliceMax = (max(max(max(handles.data(:,:,sliceNo,timeNo)))));
            sliceMin = (min(min(min(handles.data(:,:,sliceNo,timeNo)))));
        elseif handles.proj == 'sag'
            sliceMax = (max(max(max(handles.data(:,sliceNo,:,timeNo)))));
            sliceMin = (min(min(min(handles.data(:,sliceNo,:,timeNo)))));
        elseif handles.proj == 'cor'
            sliceMax = (max(max(max(handles.data(sliceNo,:,:,timeNo)))));
            sliceMin = (min(min(min(handles.data(sliceNo,:,:,timeNo)))));
        end

        if sliceMax <= sliceMin
            sliceMax = sliceMin+1;
        end
        % Update low/high number (and thus slider) values
        set( handles.lowNumber,'String',num2str(sliceMin) );
        set( handles.highNumber,'String',num2str(sliceMax) );
        set( handles.lowSlider,'Value',sliceMin );
        set( handles.highSlider,'Value',sliceMax );
        handles.highLevel = sliceMax;
        handles.lowLevel = sliceMin;
        %Change display levels
        set( handles.axes1,'CLim',[sliceMin sliceMax] );
    end
    
    % Plot image
    set( findall( handles.axes1,'type','image'), 'CData',combinedImage );
    if (datenum(handles.MATLABversion)<datenum('October 3, 2014')) %Colorbar changed after R2014b (released 3/10-2014) 
        set( handles.colorbar, 'YLim',[handles.lowLevel handles.highLevel] );
        if ~isempty(handles.colorbar)
            hCBImg = findobj(handles.colorbar,'type','image');
            set( hCBImg, 'YData',[handles.lowLevel handles.highLevel] );
        end
    else
        if ~isempty(handles.colorbar)
            set( handles.colorbar, 'Limits',[handles.lowLevel handles.highLevel] );
%             set( handles.colorbar, 'LimitsMode','auto' );
            colormap( handles.colorbar, handles.colormaps.String{get(handles.colormaps,'Value')} );
            caxis(handles.axesColorbar,[handles.lowLevel handles.highLevel]);
        end
    end
   
    % Plot contour of ROI
    if ( max( current2DROIFlat(:) ) > 0 ) && ( get(handles.hideROIcheckbox,'Value') == 0 )
        set( handles.axes1, 'NextPlot', 'add' );
        delete( findall( handles.axes1,'Tag','Contour') );
        
        for i = 1:numel(handles.ROIs)
            [~,ch1] = contour( handles.axes1,current2DROI{i}, [0.5 0.5],'LineStyle','-','LineColor',[0 0 0],'LineWidth',2);
            [~,ch2] = contour( handles.axes1,current2DROI{i}, [0.5 0.5],'LineStyle','-','LineColor',...
                handles.ROIColor{mod(i,numel(handles.ROIColor))+1},'LineWidth',1);
            set( ch1,'Tag','Contour' );
            set( ch2,'Tag','Contour' );
        end
        set( handles.axes1, 'NextPlot', 'replace' );
    else
        delete( findall( handles.axes1,'Tag','Contour') );
    end
    
end


function clearROI_Callback(hObject, eventdata, handles)
    % hObject    handle to clearROI (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    ROINumber       = get(handles.ROINumberMenu,'Value');
    noInitLegends   = numel( handles.ROILegend );

    % Drop down list
    handles.ROIs(ROINumber) = [];
    tempStr                 = cellstr( get(handles.ROINumberMenu, 'String') );
    if ( numel(handles.ROIs) == 0 )
        tempStr = {'No ROIs...'};
        set(handles.ROINumberMenu,'Enable', 'off');
    else
        tempStr = [ tempStr(1:ROINumber-1) ; tempStr(ROINumber+1:end) ];
    end
    set( handles.ROINumberMenu, 'String', tempStr );
    
    %Reset ROI statistics texts
    if ( numel(handles.ROIs) == 0 )
        set( handles.textROIName, 'String', 'ROI' );
        set( handles.textROIMean,'String','' );
        set( handles.textROIStd,'String','' );
        set( handles.textROIMin,'String','' );
        set( handles.textROIMax,'String','' );
        set( handles.textROISum,'String','' );
        set( handles.textROINoVoxels,'String','' );
    end
    
    % Delete all text legends
    for i = noInitLegends:-1:1
        try
            delete( handles.ROILegend{i} );
            handles.ROILegend(i) = [];
%             delete( handles.ROILegend{i} );
        end
    end
    
    % Write legends in text boxes
    legendStr       = get( handles.ROINumberMenu, 'String');
    for i = 1:numel(handles.ROIs)
        if i > 1
            box     = get( handles.ROILegend{i-1},'Extent' );
            yPos    = box(2)-handles.ROILegendHeight;
        else
            yPos    = handles.yMax;
        end
        handles.ROILegend{i} = text( 'Position', [1,yPos,0], 'String', legendStr{i}, 'FontSize', 12,'Color',...
            handles.ROIColor{mod(i,numel(handles.ROIColor))+1},...
            'HorizontalAlignment','Left','VerticalAlignment','bottom',...
            'Interpreter','none');
    end
    
    % Choose new element in drop down list
    if ( (ROINumber == numel( handles.ROIs ) + 1) && (ROINumber ~= 1) ) % Last element removed
        set( handles.ROINumberMenu, 'Value', numel( handles.ROIs ) );
    else 
        set( handles.ROINumberMenu, 'Value', ROINumber );
    end
    ROINumberMenu_Callback(handles.ROINumberMenu, [], handles)
    
    %Redraw image
    handles = drawImagePlusROI(handles);
    guidata(hObject,handles)
end


function ROINumberMenu_Callback(hObject, eventdata, handles)

    %If no ROI is chosen, return
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    ROINumber = get(hObject,'Value');
    ROINames  = get(hObject,'String');

    set( handles.textROIName, 'String', ROINames(ROINumber) );
    
    % Highlight chosen ROI legend
    for i = 1:numel( handles.ROILegend )
        if i == ROINumber
            set( handles.ROILegend{i}, 'FontWeight','Bold');
            set( handles.ROILegend{i}, 'FontSize',15);
        else
            set( handles.ROILegend{i}, 'FontWeight','Normal');
            set( handles.ROILegend{i}, 'FontSize',12);
            
        end
    end

    % Redra image and calculate ROI statistics
    handles   = drawImagePlusROI(handles);
    handles   = calculateROIStatistics(handles);

    guidata(hObject,handles)
end


% --- Executes on button press in addROI.
function addROI_Callback(hObject, eventdata, handles)
    % hObject    handle to addROI (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % New ROI element
    handles.ROIs{end+1} = zeros(size(handles.data,1),size(handles.data,2),size(handles.data,3));
    
    handles = addROICont(handles);
    
    guidata(hObject,handles);
end


function handles = addROICont(handles)

%     % New ROI element
%     handles.ROIs{end+1} = zeros(size(handles.data,1),size(handles.data,2),size(handles.data,3));
    
    ROINumber           = numel( handles.ROIs ); 
    set(handles.ROINumberMenu,'Enable', 'on');
    
    % Drop down list              
    tempStr             = cellstr( get(handles.ROINumberMenu, 'String') );
    if ( numel(handles.ROIs) == 1 ) && ( strcmp(tempStr(1),'No ROIs...')==1 )
        tempStr         = cell(0);
    end
    tempStr             = [ tempStr ; cellstr( get( handles.textROIName,'String' )) ];
    set( handles.ROINumberMenu, 'String', tempStr );
    set( handles.ROINumberMenu, 'Value', numel(handles.ROIs) );
    
    % Text in axes
    ROILegend   = get( handles.ROINumberMenu, 'String');
    if ROINumber > 1
        box     = get( handles.ROILegend{ROINumber-1},'Extent' );
        yPos    = box(2)-handles.ROILegendHeight;
    else
        yPos    = handles.yMax;
    end
    handles.ROILegend{ROINumber} = text( 'Position', [1,yPos,0],...
                'String', ROILegend{ROINumber},'FontSize', 12,...
                'Color',handles.ROIColor{mod(ROINumber,numel(handles.ROIColor))+1},...
                'HorizontalAlignment','Left','VerticalAlignment','bottom',...
                'Interpreter','none');
    temp                         = get(  handles.ROILegend{ROINumber}, 'Extent');
    handles.ROILegendHeight      = temp(4);
    
    % Callback for drop down list change
    ROINumberMenu_Callback(handles.ROINumberMenu, [], handles)
                             
    guidata(handles.addROI,handles)
end


% --- Executes on button press in ROIChangeName.
function ROIChangeName_Callback(hObject, eventdata, handles)
    % hObject    handle to ROIChangeName (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    ROINumber           = get(handles.ROINumberMenu,'Value');
    ROINames            = get(handles.ROINumberMenu,'String');
    ROINames(ROINumber) = strtrim(cellstr(get(handles.textROIName,'String')));
    set(handles.ROINumberMenu,'String',ROINames)

    % Text in axes
    ROILegend = get( handles.ROINumberMenu, 'String');
    set( handles.ROILegend{ROINumber}, 'String', ROILegend{ROINumber});
   
    guidata(hObject,handles)
end


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    coordinates = round( get(handles.axes1,'CurrentPoint') );
    if ( coordinates(1,1) >= handles.xMin ) && ( coordinates(1,1) <= handles.xMax ) ...
            && ( coordinates(1,2) >= handles.yMin ) && ( coordinates(1,2) <= handles.yMax )
        
        t       = str2num( get(handles.timeNumber,'String') );
        slice   = str2num( get(handles.sliceNumber,'String') );
        x       = coordinates(1,1); % horizontal in displayed image
        y       = coordinates(1,2); % vertical in displayed image
        if handles.proj == 'tra'
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(1) size( handles.data,2) size( handles.data,1) ];
            value     = handles.data(y,x,slice,t);
        elseif handles.proj == 'sag'
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(1) handles.voxelSize(3) size( handles.data,1) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(x,slice,y,t);
        elseif handles.proj == 'cor' 
            % Displayed image rotated 90 degrees compared to data
            pixelSize = [ handles.voxelSize(2) handles.voxelSize(3) size( handles.data,2) size( handles.data,3) ];
            y         = handles.yMax - y + 1;
            value     = handles.data(slice,x,y,t);       
        end

        xUnit       = ( x - pixelSize(3)/2 - 1/2 ) * pixelSize(1);
        yUnit       = ( y - pixelSize(4)/2 - 1/2 ) * pixelSize(2);
        if ( abs(value) >= 1e4 ) %|| ( value <= 1e-3 )
            string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value,'%.3e') ];
        else
            string      = ['(' num2str(x) ', ' num2str(y) ')  '  num2str(value,'%.4f') ];
        end
        stringUnit  = ['(' num2str(xUnit,'%.2f') ', ' num2str(yUnit,'%.2f') ') ' handles.dimUnits{1} ];
        set( handles.textPixelInfo,'String',string );
        set( handles.textPixelInfoUnit,'String',stringUnit );
        
        % Line profile coordinates
        if ~isempty( findall( handles.axes1,'type','line') );
            lineProfX = get( handles.lineProfileCenter,'XData');
            lineProfY = get( handles.lineProfileCenter,'YData');
        else
            lineProfX = [-1 -1];
            lineProfY = [-1 -1];
        end
        
        % If no ROI selected, do arrow pointer. If on line profile, do hand pointer.
        if ( coordinates(1,1) >= lineProfX(1) ) && ( coordinates(1,1) <= lineProfX(2) ) ...
            && ( coordinates(1,2) >= lineProfY(1) ) && ( coordinates(1,2) <= lineProfY(2) )
            set(gcf,'Pointer','hand')
            set(gcf,'Pointer','hand')
        elseif strcmp(get(handles.ROINumberMenu,'Enable'),'off')
            set(gcf,'Pointer','arrow')
        else
            set(gcf,'Pointer','crosshair')
        end
        
    else
        set( handles.textPixelInfo,'String','(Column,Row) Value' );
        set( handles.textPixelInfoUnit,'String','(Column,Row) absolute' );
        set(gcf,'Pointer','arrow')
    end

    guidata(hObject,handles)
end

% --- Executes on button press in radiobuttonSquareROI.
function radiobuttonSquareROI_Callback(hObject, eventdata, handles)
    % hObject    handle to radiobuttonSquareROI (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in radiobuttonCircleROI.
function radiobuttonCircleROI_Callback(hObject, eventdata, handles)
    % hObject    handle to radiobuttonCircleROI (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
end


% --- Executes on button press in ROIExport2WSPushbutton.
function ROIExport2WSPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to ROIExport2WSPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Screen output
    if numel(handles.ROIs) == 0
        fprintf('No ROIs to export!\n\n')
        return
    end
    
    % Screen output
    disp('Exporting all ROIs to workspace...')
    
    % All ROI names
    ROInames = get( handles.ROINumberMenu,'String');

    % Export ROIs
    for i = 1:numel( handles.ROIs )
        roi     = handles.ROIs{i};
        assignin('base',['ROI_' ROInames{i}],roi);
    end

    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles)
end


% --- Executes on button press in Scale2SliceRadioButton.
function Scale2SliceRadioButton_Callback(hObject, eventdata, handles)
    % hObject    handle to Scale2SliceRadioButton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    set( hObject,'Value',1);
    set( handles.Scale2ImageRadiobutton,'Value',0);
    
    sliceNo = str2num( get(handles.sliceNumber,'String'));
    timeNo  = str2num( get(handles.timeNumber,'String'));
    
    % Image max and min value
    if handles.proj == 'tra'
        sliceMax = (max(max(max(handles.data(:,:,sliceNo,timeNo)))));
        sliceMin = (min(min(min(handles.data(:,:,sliceNo,timeNo)))));
    elseif handles.proj == 'sag'
        sliceMax = (max(max(max(handles.data(:,sliceNo,:,timeNo)))));
        sliceMin = (min(min(min(handles.data(:,sliceNo,:,timeNo)))));
    elseif handles.proj == 'cor'
        sliceMax = (max(max(max(handles.data(sliceNo,:,:,timeNo)))));
        sliceMin = (min(min(min(handles.data(sliceNo,:,:,timeNo)))));
    end
       
    if sliceMax <= sliceMin
        step = get(handles.lowSlider,'sliderStep');
        sliceMax = sliceMin + step(1);
    end
        
    % Update low/high number (and thus slider) values
    set( handles.lowNumber,'String',num2str(sliceMin) );
    set( handles.highNumber,'String',num2str(sliceMax) );
    set( handles.lowSlider,'Value',sliceMin );
    set( handles.highSlider,'Value',sliceMax );
    highNumber_Callback( handles.highNumber,[],handles);
    lowNumber_Callback( handles.lowNumber,[],handles);

    guidata(hObject,handles)
end


% --- Executes on button press in Scale2ImageRadiobutton.
function Scale2ImageRadiobutton_Callback(hObject, eventdata, handles)
    % hObject    handle to Scale2ImageRadiobutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    set( hObject,'Value',1);
    set( handles.Scale2SliceRadioButton,'Value',0);
    
    % Image max and min value
    imageMax = ( max(handles.data(:)) );
    imageMin = ( min(handles.data(:)) );
   
    if imageMax <= imageMin
        step = get(handles.lowSlider,'sliderStep');
        imageMax = imageMin + step(1);
    end
    
    % Update low/high number (and thus slider) values
    set( handles.lowNumber,'String',num2str(imageMin) );
    set( handles.highNumber,'String',num2str(imageMax) );
    set( handles.lowSlider,'Value',imageMin );
    set( handles.highSlider,'Value',imageMax );
    handles.lowLevel = imageMin;
    handles.highLevel = imageMax;
    highNumber_Callback( handles.highNumber,[],handles);
    lowNumber_Callback( handles.lowNumber,[],handles);

    guidata(hObject,handles)
end


% --- Executes on button press in showOnlyROICheckbox.
function showOnlyROICheckbox_Callback(hObject, eventdata, handles)
    % hObject    handle to showOnlyROICheckbox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    if ( get(hObject,'Value') == 1 ) 
        set( handles.axes1,'CLim',[0 1]);
    else
        set( handles.axes1,'CLim',[handles.lowLevel handles.highLevel]);
    end
    handles = drawImagePlusROI(handles);

    guidata(hObject,handles);
end


% --- Executes on button press in ROICopyToSlicesPushbutton.
function ROICopyToSlicesPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to ROICopyToSlicesPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    %If no ROI is chosen, return
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    % Screen output
    disp('Copying current ROI to all slices...')
    
    ROINumber   = get(handles.ROINumberMenu,'Value');
    sliceNo     = str2num( get(handles.sliceNumber,'String') );

    if handles.proj == 'tra'
        currentROI = squeeze( handles.ROIs{ROINumber}(:,:,sliceNo) );
        for i = 1:size( handles.data,3)
            handles.ROIs{ROINumber}(:,:,i) = currentROI;
        end
    elseif handles.proj == 'sag'
        currentROI = squeeze( handles.ROIs{ROINumber}(:,sliceNo,:) );
        for i = 1:size( handles.data,2)
            handles.ROIs{ROINumber}(:,i,:) = currentROI;
        end
    elseif handles.proj == 'cor'
        currentROI = squeeze( handles.ROIs{ROINumber}(sliceNo,:,:) );
        for i = 1:size( handles.data,1)
            handles.ROIs{ROINumber}(i,:,:) = currentROI;
        end
    end
    
    %Re-calculate ROI statistics
    handles = calculateROIStatistics(handles);
    
    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles);
end


% --- Executes on slider movement.
function ROIWindowSlider_Callback(hObject, eventdata, handles)
    % hObject    handle to ROIWindowSlider (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles.ROIwindow   = get(handles.ROIWindowSlider,'Value');
    handles             = drawImagePlusROI(handles);
    
    guidata(hObject,handles);
end


% --- Executes on button press in invertCheckbox.
function invertCheckbox_Callback(hObject, eventdata, handles)
    % hObject    handle to invertCheckbox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
%     mycmap = colormap;
    mycmap = get(handles.figure1,'Colormap');
%     mycmap = colormap( handles.colormaps.String{get(handles.colormaps,'Value')} );
%     handles.figure1.Colormap = flipud(mycmap);
%     handles.colorbar.Colormap = flipud(mycmap);
%     mycmap = get(handles.figure1,'Colormap');
    colormap(handles.figure1,flipud(mycmap))
    colormap(handles.colorbar,flipud(mycmap))
%     colormap( handles.colorbar, handles.colormaps.String{get(handles.colormaps,'Value')} );

    guidata(hObject,handles);
end


% --- Executes on button press in ROIProfPushbutton.
function ROIProfPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to ROIProfPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    %If no ROI is chosen, return
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end
    
    % Screen output
    disp('Generating ROI profile...')
    
    legendStr = '';
    cutForMean = 4; % first and last profile values cut for mean calculation
    
    % Plot figure in existing or new window
    fh = findobj('type','figure','name','ROI profile');

    % Current frame
    timeNo      = str2num( get(handles.timeNumber,'String') );
    
    %Data in ROI for each slice
    for i = 1:numel( handles.ROIs )
        tempData    = squeeze( handles.ROIs{i} ) .* ...
            squeeze( handles.data(:,:,:,timeNo) );

        switch lower(handles.proj)
            case 'tra'
                ROIedSlice(:,i) = sum( squeeze( sum(tempData,1) ),1 ) ./ ...
                    sum( squeeze( sum(handles.ROIs{i},1) ),1 );
                tempStr         = '(in Z-dir)';
            case 'cor'
                ROIedSlice(:,i) = sum( squeeze( sum(tempData,2) ),2 ) ./ ...
                    sum( squeeze( sum(handles.ROIs{i},2) ),2 );
                tempStr         = '(in X-dir)';
            case 'sag'
                ROIedSlice(:,i) = sum( squeeze( sum(tempData,3) ),1 ) ./ ...
                    sum( squeeze( sum(handles.ROIs{i},3) ),1 );
                tempStr         = '(in Y-dir)';
        end
        % Remove NaN numbers (set to zero)
        ROIedSlice( isnan(ROIedSlice) ) = 0;

        meanProf(i) = mean(ROIedSlice(cutForMean+1:end-cutForMean,i));
        if isempty( fh )
            fh = figure('name','ROI profile','numbertitle','off','units','normalized'); 
            hold on;
            
            % Plot mean value of profile
            ah = plot( [cutForMean+1 size(ROIedSlice(:,i),1)-cutForMean],[meanProf(i) meanProf(i)], ':+','Color',[.5 .5 .5]);
            % Plot profile
            ph(i) = plot( ROIedSlice(:,i), '-s','LineWidth',1,'MarkerSize',4,'Color',[0 0 0],'MarkerFaceColor',...
                handles.ROIColor{mod(i,numel(handles.ROIColor))+1},'MarkerEdgeColor',[0 0 0]);
        else
            figure(fh);
          
            % Plot mean value of profile
            ah = plot( [cutForMean+1 size(ROIedSlice(:,i),1)-cutForMean],[meanProf(i) meanProf(i)], ':+','Color',[.5 .5 .5]);
            hold on;
            % Plot profile
            ph(i) = plot( ROIedSlice(:,i), '-s','LineWidth',1,'MarkerSize',4,'Color',[0 0 0],'MarkerFaceColor',...
                 handles.ROIColor{mod(i,numel(handles.ROIColor))+1},'MarkerEdgeColor',[0 0 0]);
        end
    end
    figure(fh); hold off;
    
    xlabel(['Slice ' tempStr],'FontSize',12), ylabel(['Intensity  (' handles.imageUnit ')'],'FontSize',12);
    title('ROI profile','FontSize',12);
    xlim([1 size(ROIedSlice,1)]);
    namesStr = get( handles.ROINumberMenu,'String');
    for i = 1:numel(meanProf)
        legendStr{i} = [ namesStr{i} '  (' num2str(meanProf(i),'%.2f') ')' ];
    end
    legend( ph, legendStr );
    
    % Set focus back to ROIViewer4D
    set(handles.figure1,'Visible','on')
    
    % Check if you should export to workspace
    if ( get( handles.ROIProfExportCheckbox,'Value' ) == 1 )
        for i = 1:numel( handles.ROIs )
            assignin('base',['ROIProfile_' handles.proj '_' namesStr{i}],ROIedSlice(:,i));
        end
    end
    
    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles);
end


% --- Executes on button press in lineProfVerticalPushbutton.
function lineProfVerticalPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfVerticalPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
 
    % Screen output
    disp('Generating vertical line profile...')
    
    sliceNo = str2num( get(handles.sliceNumber,'String') );
    timeNo  = str2num( get(handles.timeNumber,'String') );
    
    switch lower(handles.proj)
        case 'tra'
            image   = squeeze( handles.data(:,:,sliceNo,timeNo) );
            col     = round( size( handles.data,2)/2 );
            x       = [col col];
            y       = [0 size( handles.data,1)+1];
        case 'cor'
            %Display rotated 90 degrees
            image   = rot90( squeeze( handles.data(sliceNo,:,:,timeNo) ), 1);
            col     = round( size( handles.data,2)/2 );
            x       = [col col];
            y       = [0 size( handles.data,3)+1];
        case 'sag'
            %Display rotated 90 degrees
            image   = rot90( squeeze( handles.data(:,sliceNo,:,timeNo) ), 1);
            col     = round( size( handles.data,1)/2 );
            x       = [col col];
            y       = [0 size( handles.data,3)+1];
    end
  
    % Size of rectagle
    width     = get( handles.sliderLineProfWidth,'Value');
    if mod(width,2) == 0
        recPosition = col - width/2 + 1;
        width       = width - 1;
    else
        recPosition = col - width/2;
    end
    
    % Draw line in image
    delete( findall( handles.axes1,'type','line') );
    delete( findall( handles.axes1,'type','rectangle') );
    set( handles.axes1, 'NextPlot', 'add' );
    handles.lineProfile          = rectangle( 'Parent',handles.axes1,'Position',[recPosition,0.5,width,handles.yMax-0.5],...
        'LineWidth',1);
    handles.lineProfileCenter    = plot( handles.axes1, x,y,':w','LineWidth',1);
    handles.lineProfileDirection = 'Vertical';
    set( handles.axes1, 'NextPlot', 'replace' );
    
    % Check what line color to use
    if get( handles.lineProfRedRadiobutton, 'Value') == 1
        set( handles.lineProfile,'EdgeColor',[1 0 0] ); % red
    else
        set( handles.lineProfile,'EdgeColor',[0 0 1] ); % Blue
    end
    
    % Plot figure in existing or new window
    handles = plotLineProfile(handles,image,col);  
    
    % Enable width slider
    handles = lineProfSliderRecalc(handles);
    set( handles.sliderLineProfWidth,'Enable','on');
    set( handles.textLineProfWidth,'Enable', 'on' );
    set( handles.textLineProfText,'Enable', 'on' );

    % Set focus back to ROIViewer4D
    set(handles.figure1,'Visible','on')
    
    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles);
end

% --- Executes on button press in lineProfHorizontalPushbutton.
function lineProfHorizontalPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfHorizontalPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Screen output
    disp('Generating horizontal line profile...')
    
    sliceNo = str2num( get(handles.sliceNumber,'String') );
    timeNo  = str2num( get(handles.timeNumber,'String') );
    
    switch lower(handles.proj)
        case 'tra'
            image = squeeze( handles.data(:,:,sliceNo,timeNo) );
            row   = round( size( handles.data,1)/2 );
            y     = [row row];
            x     = [0 size( handles.data,2)+1];
        case 'cor'
            %Display rotated 90 degrees
            image = rot90(squeeze( handles.data(sliceNo,:,:,timeNo) ), 1);
            row   = round( size( handles.data,3)/2 );
            y     = [row row];
            x     = [0 size( handles.data,2)+1];
        case 'sag'
            %Display rotated 90 degrees
            image = rot90(squeeze( handles.data(:,sliceNo,:,timeNo) ), 1);
            row   = round( size( handles.data,3)/2 );
            y     = [row row];
            x     = [0 size( handles.data,1)+1];
    end

    % Size of rectagle
    width       = get( handles.sliderLineProfWidth,'Value');
    if mod(width,2) == 0
        recPosition = row - width/2 + 1;
        width       = width - 1;
    else
        recPosition = row - width/2;
    end
    
    % Draw line in image
    delete( findall( handles.axes1,'type','line') );
    delete( findall( handles.axes1,'type','rectangle') );
    set( handles.axes1, 'NextPlot', 'add' );
    handles.lineProfile          = rectangle( 'Parent',handles.axes1,'Position',[0.5,recPosition,handles.xMax-0.5,width],...
        'LineWidth',1);
    handles.lineProfileCenter    = plot( handles.axes1,x,y,':w','LineWidth',1);
    handles.lineProfileDirection = 'Horizontal';
    set( handles.axes1, 'NextPlot', 'replace' );
    
    % Check what line color to use
    if get( handles.lineProfRedRadiobutton, 'Value') == 1
        set( handles.lineProfile,'EdgeColor',[1 0 0] ); % red
    else
        set( handles.lineProfile,'EdgeColor',[0 0 1] ); % Blue
    end
    
    % Plot figure in existing or new window
    handles = plotLineProfile(handles,image,row); 

    % Enable width slider
    handles = lineProfSliderRecalc(handles);
    set( handles.sliderLineProfWidth,'Enable','on');
    set( handles.textLineProfWidth,'Enable', 'on' );
    set( handles.textLineProfText,'Enable', 'on' );

    % Set focus back to ROIViewer4D
    set(handles.figure1,'Visible','on')
    
    % Screen output
    fprintf('Done!\n\n')
    
    guidata(hObject,handles);
end


function handles = lineProfSliderRecalc(handles)

    switch handles.lineProfileDirection
        case 'Vertical'
            set( handles.sliderLineProfWidth,'Max',handles.xMax);
            set( handles.sliderLineProfWidth,'SliderStep',[1/(handles.xMax-1) 10/(handles.xMax-1)])
        case 'Horizontal'
            set( handles.sliderLineProfWidth,'Max',handles.yMax);
            set( handles.sliderLineProfWidth,'SliderStep',[1/(handles.yMax-1) 10/(handles.yMax-1)])
    end
    
    valueOld = get( handles.sliderLineProfWidth, 'Value' );
    
    if valueOld > get( handles.sliderLineProfWidth,'Max' )
        set( handles.sliderLineProfWidth,'Value', get( handles.sliderLineProfWidth,'Max' ) );
    end

    set( handles.textLineProfWidth,'String', [ num2str( get( handles.sliderLineProfWidth,'Value' )) ...
        ' / ' num2str(get(handles.sliderLineProfWidth,'Max')) ] );
    
end


% --- Executes on slider movement.
function sliderLineProfWidth_Callback(hObject, eventdata, handles)
    % hObject    handle to sliderLineProfWidth (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    value   = round( get( hObject,'Value') );
    set( handles.textLineProfWidth,'String', [num2str(value) ' / ' num2str(get(hObject,'Max'))] );
    
    col     = get( handles.lineProfileCenter,'XData');
    row     = get( handles.lineProfileCenter,'YData');
    handles = drawLineProfileInAxis(handles,col(1),row(1),value);
    
    guidata(hObject,handles);
end


function  handles = plotLineProfile(handles,image,index)
    
    if strcmp( handles.lineProfileDirection,'Vertical')
        temp = mean( image(:,index), 2);
        x    = [1:numel(temp)]';
    elseif strcmp( handles.lineProfileDirection,'Horizontal')
        temp = mean( image(index,:), 1);
        x    = [1:numel(temp)];
    end
    
    fh = findobj('type','figure','name','Line profile');
    if isempty( fh )
        fh = figure('name','Line profile','numbertitle','off');
        plot(temp,'-s','LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.5 0.9 0.7],...
            'MarkerEdgeColor',[0 0 0]);
        ah = findobj( fh, 'type','axes');
    else
        ah = findobj( fh, 'type','axes');
        plot( ah, temp,'-s','LineWidth',1,'MarkerSize',4,'MarkerFaceColor',[0.5 0.9 0.7],...
            'MarkerEdgeColor',[0 0 0]);
    end
    xlabel(ah,'Pixel','FontSize',12);
    ylabel(ah,['Intensity  (' handles.imageUnit ')'],'FontSize',12);
    title(ah,'Line profile','FontSize',12);
    xlim(ah,[1 numel(temp)]);
    
end


% --- Executes on button press in lineProfExportPushbutton.
function lineProfExportPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfExportPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % If no line profile, do nothing 
    if isempty( findall( handles.axes1,'type','line') )
        return
    end
    
    sliceNo = str2num( get(handles.sliceNumber,'String') );
    timeNo  = str2num( get(handles.timeNumber,'String') );
    
    switch lower(handles.proj)
        case 'tra'
            image = squeeze( handles.data(:,:,sliceNo,timeNo) );
        case 'cor'
            image = rot90( squeeze( handles.data(sliceNo,:,:,timeNo) ) );
        case 'sag'
            image = rot90( squeeze( handles.data(:,sliceNo,:,timeNo) ) );
    end
    
    switch handles.lineProfileDirection
        case 'Vertical'
            nameStr = ['lineProfile_' handles.proj '_V'];
            index   = get( handles.lineProfileCenter, 'XData');
            temp    = image(:,index(1));
        case 'Horizontal'
            nameStr = ['lineProfile_' handles.proj '_H'];
            index   = get( handles.lineProfileCenter, 'YData');
            temp    = image(index(1),:)';
    end
    
    assignin('base', nameStr, temp );
    
    guidata(hObject,handles);
end


% --- Executes on button press in lineProfClearPushbutton.
function lineProfClearPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfClearPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    delete( findall( handles.axes1,'type','line') );
    delete( findall( handles.axes1,'type','rectangle') );
    
    % Disable width slider
    set( handles.sliderLineProfWidth,'Enable','off');
    set( handles.textLineProfWidth,'Enable', 'off' );
    set( handles.textLineProfText,'Enable', 'off' );
    set( handles.textLineProfWidth,'String', [num2str( get(handles.sliderLineProfWidth,'Value') )] );
    
    guidata(hObject,handles);
end


% --- Executes on button press in lineProfRedRadiobutton.
function lineProfRedRadiobutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfRedRadiobutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    set( hObject,'Value',1);
    set( handles.lineProfBlueRadiobutton,'Value',0);
    
    % If found line profile change color
    if ~isempty( findall( handles.axes1,'type','rectangle') )
        set( handles.lineProfile,'EdgeColor',[1 0 0])
    end
    
    guidata(hObject,handles);
end


% --- Executes on button press in lineProfRedRadiobutton.
function lineProfBlueRadiobutton_Callback(hObject, eventdata, handles)
    % hObject    handle to lineProfRedRadiobutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    set( hObject,'Value',1);
    set( handles.lineProfRedRadiobutton,'Value',0);
    
    % If found line profile change color
    if ~isempty( findall( handles.axes1,'type','rectangle') )
        set( handles.lineProfile,'EdgeColor',[0 0 1])
    end

    guidata(hObject,handles);
end

% --- Executes on button press in HistPushbutton.
function HistPushbutton_Callback(hObject, eventdata, handles)
    % hObject    handle to HistPushbutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    %If no ROI is chosen, return
    if strcmp(get(handles.ROINumberMenu,'Enable'),'off')
        return
    end

    % Screen output
    disp('Calculating histogram in current ROI...')
    
    ROINumber   = get( handles.ROINumberMenu,'Value' );
    sliceNo     = str2num( get(handles.sliceNumber,'String') );
    timeNo      = str2num( get(handles.timeNumber,'String') );
    ROINames    = get( handles.ROINumberMenu,'String' );
    noBins      = round( get( handles.sliderNoBins, 'Value') );
    
    % ROIed image
    ROIindex    = find( handles.ROIs{ROINumber} > 0 );
    try 
        if ~isempty( ROIindex )
            ROIedTemp   = handles.data(:,:,:,timeNo);
            ROIedImage  = ROIedTemp(ROIindex);
        else
            ROIedImage  = 0;
        end        
    catch
        ROIedImage  = handles.data;
    end
     
    calculateHist(ROIedImage(:),noBins,ROINames(ROINumber),handles); 
    % Screen output
    fprintf('Done!\n\n')
end

function [n,x] = calculateHist(data,noBins,ROIname,handles)

    % Histogram
    [n,x]       = hist(data,noBins); 
    
    % Plot histogram
    fh = findobj('type','figure','name','ROI histogram');
    if isempty( fh )
        fh = figure('name','ROI histogram','numbertitle','off');
        bar(x,n);     
        ah = findobj( fh, 'type','axes');
    else
        ah = findobj( fh, 'type','axes','tag','');
        bar(ah,x,n);
    end
    legend(ah,ROIname);
    xlabel(ah,['Intensity  (' handles.imageUnit ')'],'FontSize',12);
    ylabel(ah,'Number of voxels','FontSize',12);
    title(ah,'ROI histogram','FontSize',12);
    xlim(ah,[min(x)*.95 max(x)*1.05]);
    ylim(ah,[0 max(n)*1.05]);
    
end


% --- Executes on slider movement.
function sliderNoBins_Callback(hObject, eventdata, handles)
    % hObject    handle to sliderNoBins (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    noBins = round( get( hObject,'Value') );
    set( handles.noBins, 'String', noBins );
end


% --- Executes during object creation, after setting all properties.
function sliderNoBins_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to sliderNoBins (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    
    % Max and min
    set(hObject,'Min',1);
    set(hObject,'Max',100);
   
    % Set slider values
    set(hObject,'SliderStep',[1/99 10/99]);
    
end


% --- Executes on button press in pushbuttonExportFig.
function pushbuttonExportFig_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbuttonExportFig (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    [filename, filepath] = uiputfile('figure.emf', 'Save As');
  
    if  strcmp( class(filename),'char')

        newFig = figure('Position', [100 100 629 630]);
        newA = copyobj(handles.axes1,newFig);
        
        cmapNames = get(handles.colormaps,'String');
        colormap( newA, cmapNames{get(handles.colormaps,'Value')} );
        
        saveas(newFig,fullfile(filepath,filename))
        close(newFig);
    end
end

function updateMinMaxText(handles)

    % Max and min of image
    if ( abs(max(handles.data(:))) >= 1e4 )
        maxImStr = num2str( max(handles.data(:)),'%.3e');
    else
        maxImStr = num2str( max(handles.data(:)),'%.4f' );
    end
    if ( abs(min(handles.data(:))) >= 1e4 )
        minImStr = num2str( min(handles.data(:)),'%.3e' );
    else
        minImStr = num2str( min(handles.data(:)),'%.4f' );
    end

    % Chosen slice/frame min max
    sliceNo = round(get(handles.sliceSlider,'Value'));
    frameNo = round(get(handles.timeSlider,'Value'));
    switch handles.proj
        case 'tra'
            slice  = squeeze(handles.data(:,:,sliceNo,frameNo));
        case 'sag'
            slice  = squeeze(handles.data(:,sliceNo,:,frameNo));
        case 'cor'
            slice  = squeeze(handles.data(sliceNo,:,:,frameNo));
    end
    if ( abs(min(slice(:))) >= 1e4 )
        minSlStr = num2str( min(slice(:)),'%.3e' );
    else
        minSlStr = num2str( min(slice(:)),'%.4f' );
    end
    if ( abs(max(slice(:))) >= 1e4 )
        maxSlStr = num2str( max(slice(:)),'%.3e' );
    else
        maxSlStr = num2str( max(slice(:)),'%.4f' );
    end
    set(handles.textImageMin,'String',minImStr);
    set(handles.textImageMax,'String',maxImStr);
    set(handles.textSliceMin,'String',minSlStr);
    set(handles.textSliceMax,'String',maxSlStr);
    
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  structure with the following fields (see FIGURE)
    %	VerticalScrollCount: signed integer indicating direction and number of clicks
    %	VerticalScrollAmount: number of lines scrolled for each click
    % handles    structure with handles and user data (see GUIDATA)
    
    % If inside the image
    coordinates = round( get(handles.axes1,'CurrentPoint') );
    if ( coordinates(1,1) >= handles.xMin ) && ( coordinates(1,1) <= handles.xMax ) ...
            && ( coordinates(1,2) >= handles.yMin ) && ( coordinates(1,2) <= handles.yMax )
    
        direction    = eventdata.VerticalScrollCount;
        scrollAmount = eventdata.VerticalScrollAmount/3;
    
        % Scroll slices
        if handles.sliceSlider.Max>1
            currentVal   = get(handles.sliceSlider,'Value');
            newValue     = currentVal-(scrollAmount*direction);
            % Inside alowed slice interval
            if newValue >= 1 && newValue <= get(handles.sliceSlider,'Max')
                set(handles.sliceSlider,'Value', currentVal-(scrollAmount*direction) );
                sliceSlider_Callback( handles.sliceSlider, [], handles )
                guidata(hObject,handles)
            end
        % Scroll frames
        elseif handles.timeSlider.Max>1
            currentVal   = get(handles.timeSlider,'Value');
            newValue     = currentVal-(scrollAmount*direction);
            % Inside alowed frame interval
            if newValue >= 1 && newValue <= get(handles.timeSlider,'Max')
                set(handles.timeSlider,'Value', currentVal-(scrollAmount*direction) );
                timeSlider_Callback( handles.timeSlider, [], handles )
                guidata(hObject,handles)
            end
        end
    end
    
end
