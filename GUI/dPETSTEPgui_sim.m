function varargout = dPETSTEPgui_sim(varargin)
    % TEST MATLAB code for test.fig
    %      TEST, by itself, creates a new TEST or raises the existing
    %      singleton*.
    %
    %      H = TEST returns the handle to a new TEST or the handle to
    %      the existing singleton*.
    %
    %      TEST('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in TEST.M with the given input arguments.
    %
    %      TEST('Property','Value',...) creates a new TEST or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before test_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to test_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help test

    % Last Modified by GUIDE v2.5 10-Aug-2017 15:35:44

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @dPETSTEPgui_sim_OpeningFcn, ...
        'gui_OutputFcn',  @dPETSTEPgui_sim_OutputFcn, ...
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
end

% --- Executes just before test is made visible.
function dPETSTEPgui_sim_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user dataInputName (see GUIDATA)
    % varargin   command line arguments to test (see VARARGIN)
    
    %% Choose default command line output for test
    handles.output = hObject;
    
    %% Tab order
    buttons = {'buttonRun','buttonReset','buttonLoad','buttonSave'};
    for i = 1:numel(buttons)
        eval(['uistack(handles.' buttons{i} ',''bottom'')'])
    end
    panel = {'uipanel6','uipanel5','uipanel4','uipanel3','uipanel2','uipanel1'};
    for i = 1:numel(panel)
        eval(['uistack(handles.' panel{i} ',''bottom'')'])
        if strcmp(panel{i},'uipanel1')
            panelGroups = {'frameSource','dataSource'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
        elseif strcmp(panel{i},'uipanel2')
            uistack(handles.countSens,'bottom')
            uistack(handles.maxRingDiff,'bottom')
            uistack(handles.blurT,'bottom')
            uistack(handles.numAngBin,'bottom')
            uistack(handles.radialBinSizeSource,'bottom')  
        elseif strcmp(panel{i},'uipanel4')
            panelGroups = {'addBioVariability','Interp','KineticModel','inputFunction'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
        elseif strcmp(panel{i},'uipanel5')
            eval(['uistack(handles.' panel{i} ',''bottom'')'])
            panelGroups = {'Postfiltering','ReconMethod'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
            uistack(handles.fovSize,'bottom')
            uistack(handles.simSize,'bottom')
        end
    end

    %% Initial toggle buttons
    handles.zfilterFromWS.Value = 1;

    %% Update handles structure
    guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = dPETSTEPgui_sim_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user dataInputName (see GUIDATA)
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inputPIM_Callback(hObject, eventdata, handles)

    %% Unpush dynamic image button
    handles.inputDYN.Value = 0;
    
    %% Enable kinetic modelling panel
    objs   = allchild(handles.uipanel4);
    styles = {'text','edit','radiobutton','checkbox'};
    for i=1:numel(styles)
        obj = findobj(objs,'style',styles{i});
        for j=1:numel(obj)
            obj(j).Enable = 'on';
        end
    end
    % Keep some fields disabled (maybe)
    if handles.SumExp.Value==0
        handles.numExp.Enable     = 'off'; 
        handles.textNumExp.Enable = 'off';
    end
    if handles.addBioVar.Value==0
        handles.varScale.Enable     = 'off'; 
        handles.textVarScale.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end


function inputDYN_Callback(hObject, eventdata, handles)

    %% Unpush dynamic image button
    handles.inputPIM.Value = 0;
    
    %% Disable kinetic modelling panel
    objs   = allchild(handles.uipanel4);
    styles = {'text','edit','radiobutton','checkbox'};
    for i=1:numel(styles)
        obj = findobj(objs,'style',styles{i});
        for j=1:numel(obj)
            obj(j).Enable = 'off';
        end
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataFromWS_Callback(hObject, eventdata, handles)
    
    %% Reset color
    set(handles.dataInputName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function dataFromFile_Callback(hObject, eventdata, handles)

    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat'},'Pick a file');
    
    %% Reset color
    set(handles.dataInputName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.dataInputName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function frameFromWS_Callback(hObject, eventdata, handles)

    %% Reset color
    set(handles.frameName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function frameFromFile_Callback(hObject, eventdata, handles)

    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat;*.txt;*.csv*'},'Pick a file');
    
    %% Reset color
    set(handles.frameName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.frameName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function dataInputName_Callback(hObject, eventdata, handles)

    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
end


function frameName_Callback(hObject, eventdata, handles)

    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fromRadBin_Callback(hObject, eventdata, handles)

    %% Enable number of angular bins
    handles.numRadBin.Enable          = 'on';
    handles.textNumRadBin.Enable      = 'on';
    handles.textNumRadBinExtra.Enable = 'on';
    
    %% Disable ring diameter
    handles.ringDiameter.Enable         = 'off';
    handles.textRingDiameter.Enable     = 'off';
    handles.textRingDiameterUnit.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end


function fromScannerDiam_Callback(hObject, eventdata, handles)

    %% Disable number of angular bins
    handles.numRadBin.Enable          = 'off';
    handles.textNumRadBin.Enable      = 'off';
    handles.textNumRadBinExtra.Enable = 'off';
    
    %% Enable ring diameter
    handles.ringDiameter.Enable         = 'on';
    handles.textRingDiameter.Enable     = 'on';
    handles.textRingDiameterUnit.Enable = 'on';
    
    %% Update handles structure
    guidata(hObject, handles);
end


function numRadBin_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function ringDiameter_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function numAngBin_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function blurT_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfNonNeg(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function maxRingDiff_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function countSens_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function noREP_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function scaleFactor_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function SF_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfNonNeg(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function RF_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfNonNeg(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function halflife_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function noDecay_Callback(hObject, eventdata, handles)

    %% Disable halflife box or not
    if hObject.Value==1 %Box ticked - no decay
        handles.halflife.Enable = 'off';
        handles.textHalflife.Enable = 'off';
        handles.textHalflifeUnit.Enable = 'off';
    else %Box unticked - decay
        handles.halflife.Enable = 'on';
        handles.textHalflife.Enable = 'on';
        handles.textHalflifeUnit.Enable = 'on';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ifFromWS_Callback(hObject, eventdata, handles)
    %% Reset color
    set(handles.inputFuncName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end


function ifFromFile_Callback(hObject, eventdata, handles)
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat;*.txt;*.csv*'},'Pick a file');
    
    %% Reset color
    set(handles.inputFuncName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.inputFuncName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
    
end


function inputFuncName_Callback(hObject, eventdata, handles)

    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function CifScaleFactor_Callback(hObject, eventdata, handles)
    
    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function OneTissue_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable = 'off';
    handles.textNumExp.Enable = 'off';
       
    %% Update handles structure
    guidata(hObject, handles);
end


function TwoTissue_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable = 'off';
    handles.textNumExp.Enable = 'off';
       
    %% Update handles structure
    guidata(hObject, handles);
end


function FRTM_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable = 'off';
    handles.textNumExp.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end


function SRTM_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable = 'off';
    handles.textNumExp.Enable = 'off';
       
    %% Update handles structure
    guidata(hObject, handles);
end


function SumExp_Callback(hObject, eventdata, handles)

    %% Enable number of exponentials box
    handles.numExp.Enable     = 'on';
    handles.textNumExp.Enable = 'on';
    
    %% Update handles structure
    guidata(hObject, handles);
end


function numExp_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function linear_Callback(hObject, eventdata, handles)

    %% Update handles structure
    guidata(hObject, handles);
end


function nearest_Callback(hObject, eventdata, handles)

    %% Update handles structure
    guidata(hObject, handles);
end


function pchip_Callback(hObject, eventdata, handles)

    %% Update handles structure
    guidata(hObject, handles);
end


function spline_Callback(hObject, eventdata, handles)

    %% Update handles structure
    guidata(hObject, handles);
end


function timeStep_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

function noBioVar_Callback(hObject, eventdata, handles)
    %% Disable box or not
    handles.varScale.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function addBioVar_Callback(hObject, eventdata, handles)
    %% Disable box or not
    handles.varScale.Enable = 'on';
    
    %% Update handles structure
    guidata(hObject, handles);
end


function varScale_Callback(hObject, eventdata, handles)
    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simSize_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

function fovSize_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

function doFBP_Callback(hObject, eventdata, handles)

    %%  Enable postfiltering
    if hObject.Value==1 || handles.doOSEM.Value==1 || handles.doOSEMPSF.Value==1
        handles.checkPostXY.Enable = 'on';
        handles.checkPostZ.Enable  = 'on';
        if handles.checkPostXY.Value == 1;
            handles.postfilterXY.Enable          = 'on';
            handles.textPostfilterXY.Enable      = 'on';
            handles.textPostfilterXYUnit.Enable  = 'on';
        end
        if handles.checkPostZ.Value == 1;
            handles.postfilterZ.Enable       = 'on';
            handles.textPostfilterZ.Enable   = 'on';
            handles.zfilterFromWS.Enable     = 'on';
            handles.zfilterFromFile.Enable   = 'on';
        end
    else
        handles.checkPostXY.Enable = 'off';
        handles.checkPostZ.Enable  = 'off';
        handles.postfilterXY.Enable          = 'off';
        handles.textPostfilterXY.Enable      = 'off';
        handles.textPostfilterXYUnit.Enable  = 'off';
        handles.postfilterZ.Enable       = 'off';
        handles.textPostfilterZ.Enable   = 'off';
        handles.zfilterFromWS.Enable     = 'off';
        handles.zfilterFromFile.Enable   = 'off';
    end
    
    %% Enable output
    if hObject.Value==1
        handles.FBP.Enable     = 'on';
        handles.FBP.Value      = 1;
        handles.FBPName.Enable = 'on';
    else
        handles.FBP.Enable     = 'off';
        handles.FBP.Value      = 0;
        handles.FBPName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end


function doOSEM_Callback(hObject, eventdata, handles)

    %% Disable iteration/subset box or not
    if hObject.Value==1 || handles.doOSEMPSF.Value==1 %Box ticked
        handles.iterNUM.Enable     = 'on';
        handles.subNUM.Enable      = 'on';
        handles.textIterNUM.Enable = 'on';
        handles.textSubNUM.Enable  = 'on';
        
        handles.checkPostXY.Enable = 'on';
        handles.checkPostZ.Enable  = 'on';
        if handles.checkPostXY.Value == 1;
            handles.postfilterXY.Enable          = 'on';
            handles.textPostfilterXY.Enable      = 'on';
            handles.textPostfilterXYUnit.Enable  = 'on';
        else
            handles.postfilterXY.Enable          = 'off';
            handles.textPostfilterXY.Enable      = 'off';
            handles.textPostfilterXYUnit.Enable  = 'off';
        end
        if handles.checkPostZ.Value == 1;
            handles.postfilterZ.Enable       = 'on';
            handles.textPostfilterZ.Enable   = 'on';
            handles.zfilterFromWS.Enable     = 'on';
            handles.zfilterFromFile.Enable   = 'on';
        else
            handles.postfilterZ.Enable       = 'off';
            handles.textPostfilterZ.Enable   = 'off';
            handles.zfilterFromWS.Enable     = 'off';
            handles.zfilterFromFile.Enable   = 'off';
        end
    else %Box unticked
        handles.iterNUM.Enable     = 'off';
        handles.subNUM.Enable      = 'off';
        handles.textIterNUM.Enable = 'off';
        handles.textSubNUM.Enable  = 'off';
        if handles.doFBP.Value ==0
            handles.checkPostXY.Enable = 'off';
            handles.checkPostZ.Enable  = 'off';
            handles.postfilterXY.Enable          = 'off';
            handles.textPostfilterXY.Enable      = 'off';
            handles.textPostfilterXYUnit.Enable  = 'off';
            handles.postfilterZ.Enable       = 'off';
            handles.textPostfilterZ.Enable   = 'off';
            handles.zfilterFromWS.Enable     = 'off';
            handles.zfilterFromFile.Enable   = 'off';
        else %only FBP
            handles.checkPostXY.Enable = 'on';
            handles.checkPostZ.Enable  = 'on';
            if handles.checkPostXY.Value == 1;
                handles.postfilterXY.Enable          = 'on';
                handles.textPostfilterXY.Enable      = 'on';
                handles.textPostfilterXYUnit.Enable  = 'on';
            else
                handles.postfilterXY.Enable          = 'off';
                handles.textPostfilterXY.Enable      = 'off';
                handles.textPostfilterXYUnit.Enable  = 'off';
            end
            if handles.checkPostZ.Value == 1;
                handles.postfilterZ.Enable       = 'on';
                handles.textPostfilterZ.Enable   = 'on';
                handles.zfilterFromWS.Enable     = 'on';
                handles.zfilterFromFile.Enable   = 'on';
            else
                handles.postfilterZ.Enable       = 'off';
                handles.textPostfilterZ.Enable   = 'off';
                handles.zfilterFromWS.Enable     = 'off';
                handles.zfilterFromFile.Enable   = 'off';
            end
        end
    end
    
    %% Enable output
    if hObject.Value==1
        handles.OSEM.Enable     = 'on';
        handles.OSEM.Value      = 1;
        handles.OSEMName.Enable = 'on';
    else
        handles.OSEM.Enable     = 'off';
        handles.OSEM.Value      = 0;
        handles.OSEMName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);

end


function doOSEMPSF_Callback(hObject, eventdata, handles)

    %% Disable iteration/subset box or not
    if hObject.Value==1 || handles.doOSEM.Value==1 %Box ticked
        handles.iterNUM.Enable     = 'on';
        handles.subNUM.Enable      = 'on';
        handles.textIterNUM.Enable = 'on';
        handles.textSubNUM.Enable  = 'on';
        
        handles.checkPostXY.Enable = 'on';
        handles.checkPostZ.Enable  = 'on';
        if handles.checkPostXY.Value == 1;
            handles.postfilterXY.Enable          = 'on';
            handles.textPostfilterXY.Enable      = 'on';
            handles.textPostfilterXYUnit.Enable  = 'on';
        else
            handles.postfilterXY.Enable          = 'off';
            handles.textPostfilterXY.Enable      = 'off';
            handles.textPostfilterXYUnit.Enable  = 'off';
        end
        if handles.checkPostZ.Value == 1;
            handles.postfilterZ.Enable       = 'on';
            handles.textPostfilterZ.Enable   = 'on';
            handles.zfilterFromWS.Enable     = 'on';
            handles.zfilterFromFile.Enable   = 'on';
        else
            handles.postfilterZ.Enable       = 'off';
            handles.textPostfilterZ.Enable   = 'off';
            handles.zfilterFromWS.Enable     = 'off';
            handles.zfilterFromFile.Enable   = 'off';
        end
    else %Box unticked
        handles.iterNUM.Enable = 'off';
        handles.subNUM.Enable  = 'off';
        handles.textIterNUM.Enable = 'off';
        handles.textSubNUM.Enable  = 'off';
        if handles.doFBP.Value == 0
            handles.checkPostXY.Enable = 'off';
            handles.checkPostZ.Enable  = 'off';
            handles.postfilterXY.Enable          = 'off';
            handles.textPostfilterXY.Enable      = 'off';
            handles.textPostfilterXYUnit.Enable  = 'off';
            handles.postfilterZ.Enable       = 'off';
            handles.textPostfilterZ.Enable   = 'off';
            handles.zfilterFromWS.Enable     = 'off';
            handles.zfilterFromFile.Enable   = 'off';
        else %only FBP
            handles.checkPostXY.Enable = 'on';
            handles.checkPostZ.Enable  = 'on';
            if handles.checkPostXY.Value == 1;
                handles.postfilterXY.Enable          = 'on';
                handles.textPostfilterXY.Enable      = 'on';
                handles.textPostfilterXYUnit.Enable  = 'on';
            else
                handles.postfilterXY.Enable          = 'off';
                handles.textPostfilterXY.Enable      = 'off';
                handles.textPostfilterXYUnit.Enable  = 'off';
            end
            if handles.checkPostZ.Value == 1;
                handles.postfilterZ.Enable       = 'on';
                handles.textPostfilterZ.Enable   = 'on';
                handles.zfilterFromWS.Enable     = 'on';
                handles.zfilterFromFile.Enable   = 'on';
            else
                handles.postfilterZ.Enable       = 'off';
                handles.textPostfilterZ.Enable   = 'off';
                handles.zfilterFromWS.Enable     = 'off';
                handles.zfilterFromFile.Enable   = 'off';
            end
        end
    end
    if hObject.Value==1 %Box ticked
        handles.psf.Enable         = 'on';
        handles.textPsf.Enable     = 'on';
        handles.textPsfUnit.Enable = 'on';
    else %Box unticked
        handles.psf.Enable         = 'off';
        handles.textPsf.Enable     = 'off';
        handles.textPsfUnit.Enable = 'off';
    end

    %% Enable output
    if hObject.Value==1
        handles.OSEMpsf.Enable     = 'on';
        handles.OSEMpsf.Value      = 1;
        handles.OSEMpsfName.Enable = 'on';
    else
        handles.OSEMpsf.Enable     = 'off';
        handles.OSEMpsf.Value      = 0;
        handles.OSEMpsfName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function iterNUM_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function subNUM_Callback(hObject, eventdata, handles)
 
    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Round value to integer
            val            = round( hObject.Value );
            hObject.Value  = val;
            hObject.String = num2str(val,'%.0f');
            
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end


function psf_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfNonNeg(hObject);
        end
        %% Update handles structure
        guidata(hObject, handles);
    end  
end

function checkPostXY_Callback(hObject, eventdata, handles)
    %% Enable postfiltering maybe
    if hObject.Value == 1;
        handles.postfilterXY.Enable          = 'on';
        handles.textPostfilterXY.Enable      = 'on';
        handles.textPostfilterXYUnit.Enable  = 'on';
    else
        handles.postfilterXY.Enable          = 'off';
        handles.textPostfilterXY.Enable      = 'off';
        handles.textPostfilterXYUnit.Enable  = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function checkPostZ_Callback(hObject, eventdata, handles)

    %% Enable postfiltering maybe
    if hObject.Value == 1;
        handles.postfilterZ.Enable       = 'on';
        handles.textPostfilterZ.Enable   = 'on';
        handles.zfilterFromWS.Enable     = 'on';
        handles.zfilterFromFile.Enable   = 'on';
    else
        handles.postfilterZ.Enable       = 'off';
        handles.textPostfilterZ.Enable   = 'off';
        handles.zfilterFromWS.Enable     = 'off';
        handles.zfilterFromFile.Enable   = 'off';
    end
        
    %% Update handles structure
    guidata(hObject, handles);
end

function postfilterXY_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Update value
        set(hObject,'Value', str2num(get(hObject,'String') ) );
        
        %% Check if numeric value
        [hObject,error] = checkIfNumeric(hObject);
        
        if ~error
            %% Check if >0
            hObject = checkIfAboveZero(hObject);
        end
        
        %% Update handles structure
        guidata(hObject, handles);
    end    
end

function postfilterZ_Callback(hObject, eventdata, handles)

    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function zfilterFromWS_Callback(hObject, eventdata, handles)
    %% Reset color
    set(handles.postfilterZ,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end
 
function zfilterFromFile_Callback(hObject, eventdata, handles)
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat;*.txt;*.csv*'},'Pick a file');
    
    %% Reset color
    set(handles.postfilterZ,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.postfilterZ.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_Callback(hObject, eventdata, handles)
    %% Enable output
    if hObject.Value==1
        handles.dataName.Enable = 'on';
    else
        handles.dataName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function dataName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function counts_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function countsName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function countsNoise_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function countsNoiseName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function nFWprompts_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function nFWpromptsName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function FWtrues_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function FWtruesName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function FWscatters_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function FWscattersName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function FWrandoms_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function FWrandomsName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function wcc_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function wccName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function FBP_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function FBPName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function OSEM_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function OSEMName_Callback(hObject, eventdata, handles)

    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end


function OSEMpsf_Callback(hObject, eventdata, handles)
    %% Enable output
    str = [hObject.Tag 'Name'];
    if hObject.Value==1
        handles.(str).Enable = 'on';
    else
        handles.(str).Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function OSEMpsfName_Callback(hObject, eventdata, handles)
    %% Color not faulty
    set(hObject,'BackgroundColor',[1 1 1]);
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
    % Key: name of the key that was pressed, in lower case
    % Character: character interpretation of the key(s) that was pressed
    % Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
    % handles    structure with handles and user dataInputName (see GUIDATA)
    %% If enter was pressed
    if strcmp(eventdata.Key,'return')
        ch = gco; %currently selected object
        switch get(ch,'Style')
            case 'radiobutton'
                set(ch,'Value',1)
            case 'edit'
                f = get(ch,'Callback');
                f(ch,[]);
            case 'pushbutton'
                f = get(ch,'Callback');
                f(ch,[]);
            case 'checkbox'
                value   = get(ch,'Value');
                newVale = abs(value-1); %0-->1 and 1-->0
                set(ch,'Value',newVale);    
        end
        
        % Press tab to go to next
        robot = java.awt.Robot;
        robot.keyPress (java.awt.event.KeyEvent.VK_TAB); % press "tab" key
        
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function buttonSave_Callback(hObject, eventdata, handles)

    %% Save settings
    gui_saveSettings(handles)
end


function buttonLoad_Callback(hObject, eventdata, handles)

    %% Load settings
    handles = gui_loadSettings(handles,'sim');
    
    %% Update handles structure
    guidata(hObject, handles); 
end


function buttonReset_Callback(hObject, eventdata, handles)
    
    %% Reset all values
    % Reset all backgrounds to white
    editBoxes = findobj(dPETSTEPgui_sim,'style','edit');
    for i=1:numel(editBoxes)
        handles.(editBoxes(i).Tag).BackgroundColor = [1 1 1];
    end
    
    % Radio buttons
    set(handles.dataFromWS,'Value',1);
    set(handles.frameFromWS,'Value',1);
    set(handles.fromRadBin,'Value',1);
        fromRadBin_Callback(handles.fromRadBin,[], handles)
    set(handles.ifFromWS,'Value',1);
    set(handles.OneTissue,'Value',1);
        OneTissue_Callback(handles.OneTissue,[], handles)
    set(handles.linear,'Value',1);
        linear_Callback(handles.linear,[], handles)
        
    % Checkboxes
    checkBoxes = findobj(allchild(handles.uipanel6),'style','checkbox');
    for i=1:numel(checkBoxes)
        set(handles.(checkBoxes(i).Tag),'Value',1);
    end
    set(handles.noDecay,'Value',1);
        noDecay_Callback(handles.noDecay,[], handles)
    set(handles.noBioVar,'Value',1);
        noBioVar_Callback(handles.noBioVar,[], handles)
    set(handles.doFBP,'Value',1);
        doFBP_Callback(handles.doFBP,[], handles)
    set(handles.doOSEM,'Value',0);
        doOSEM_Callback(handles.doOSEM,[], handles)
    set(handles.doOSEMPSF,'Value',0);
        doOSEMPSF_Callback(handles.doOSEMPSF,[], handles)
    
    % Toggle button
    set(handles.inputPIM,'Value',1);
        inputPIM_Callback(handles.inputPIM,[], handles)
    
    % Edit boxes
    editBoxes = findobj(dPETSTEPgui_sim,'style','edit');
    for i=1:numel(editBoxes)
        if strcmp( handles.(editBoxes(i).Tag).Enable,'on')
%             handles.(editBoxes(i).Tag).String = '0';
%             handles.(editBoxes(i).Tag).Value  = 0;
            handles.(editBoxes(i).Tag).String = '';
            handles.(editBoxes(i).Tag).Value  = [];
        else
            handles.(editBoxes(i).Tag).String = '';
            handles.(editBoxes(i).Tag).Value  = [];
        end
    end
    editBoxes = findobj(allchild(handles.uipanel6),'style','edit');
    for i=1:numel(editBoxes)
        set(handles.(editBoxes(i).Tag),'String',editBoxes(i).Tag(1:end-4));
        set(handles.(editBoxes(i).Tag),'Value',[]);
    end
    set(handles.dataInputName,'String','data');
    % Special cases
    set(handles.timeStep,'Value',1);       set(handles.timeStep,'String',1);
    set(handles.CifScaleFactor,'Value',1); set(handles.CifScaleFactor,'String',1);
    set(handles.noREP,'Value',1);          set(handles.noREP,'String',1);
    set(handles.dataInputName,'Value',[]); set(handles.dataInputName,'String','data');
    set(handles.frameName,'Value',[]);     set(handles.frameName,'String','frame');
    set(handles.inputFuncName,'Value',[]); set(handles.inputFuncName,'String','Cif');
    set(handles.postfilterZ,'Value',[]);   set(handles.postfilterZ,'String','postfilterZ');

    %% Update handles structure
    guidata(hObject, handles);
end


function buttonRun_Callback(hObject, eventdata, handles)

    %% Current pointer
    oldpointer = get(handles.figure1,'pointer');
    if strcmp(oldpointer,'watch')
        oldpointer = 'arrow';
    end
    
    %% If pressing ctrl+c (abort)
    funcAbort = onCleanup(@() myCleanupFun(handles,oldpointer));

    %% Check that variable name input exist on WS or as file
    error      = 0;
    WS         = evalin('base','whos');
    varName    = {handles.dataInputName   handles.frameName   handles.inputFuncName handles.postfilterZ};
    toggleName = {handles.dataFromWS handles.frameFromWS handles.ifFromWS      handles.zfilterFromWS}; 
    for i=1:numel(varName)
        
        if ( toggleName{i}.Value==1 || isempty(toggleName{i}) ) && ( strcmp(varName{i}.Enable,'on') ) %selected read from WS
            
            if ~ismember(varName{i}.String,{WS(:).name})
                fprintf('Variable "%s" is not an existing WS variable.\n',varName{i}.String)
                set(varName{i},'BackgroundColor',[1 0.4 0.4]);
                error = 1;
            else
                tmp   = evalin('base',varName{i}.String);
                if isempty(tmp)
                    fprintf('Variable "%s" is empty.\n',varName{i}.String);
                    set(varName{i},'BackgroundColor',[1 0.4 0.4]);
                else
                    set(varName{i},'BackgroundColor',[1 1 1]);
                end
            end
           
        elseif strcmp(varName{i}.Enable,'on') %selected read from file
            if exist(varName{i}.String, 'file')~=2
                fprintf('File "%s" does not exist.\n',varName{i}.String);
                set(varName{i},'BackgroundColor',[1 0.4 0.4]);
                error = 1;
            else
                [~,~,ext] = fileparts(varName{i}.String);
                if strcmp(ext,'.mat') %MAT-file
                    tmp = cell2mat(struct2cell( load(varName{i}.String) ));
                else %Text file
                    tmp = load(varName{i}.String);
                end
                if ~isempty(tmp)
                    set(varName{i},'BackgroundColor',[1 1 1]);
                else
                    fprintf('File "%s" is empty.\n',varName{i}.String);
                    set(varName{i},'BackgroundColor',[1 0.4 0.4]);
                end
            end
        end

    end
    
    %% Check that all user input is ok
    boxes = findobj(dPETSTEPgui_sim,'style','edit');
    for i=1:numel(boxes)
        str   = eval(['handles.' boxes(i).Tag '.String']);
        color = eval(['handles.' boxes(i).Tag '.BackgroundColor']);
        on    = eval(['handles.' boxes(i).Tag '.Enable']);
        if (~isequal(color,[1 1 1]) || isempty(str) ) && strcmp(on,'on')
            obj =  eval(['handles.' boxes(i).Tag]);
            fprintf('Box "%s" is faulty.\n',get(obj,'Tag'));
            set(obj,'BackgroundColor',[1 0.4 0.4]);
            error = 1;
        end
    end
       
    %% Break if error
    if error; return; end
    
    %% Update handles structure
    guidata(hObject, handles);
      
    %% Indicate running with pointer
    handles.buttonRun.Value = 1;
    set(handles.figure1,'pointer','watch');  
    handles.running.String  = 'Running...';
    drawnow;
    
    %% Read data from WS or file
    if handles.dataFromWS.Value==1 %selected read from WS
        data  = evalin('base',handles.dataInputName.String);
    else
        data  = cell2mat(struct2cell( load(handles.dataInputName.String) ));
    end

    %% Put all settings into variable
    simSet   = populateSimSet(handles);
    types    = {data.type};
    for i=1:numel(types)
        if ~isempty(types{i})
            switch lower(types{i})
                case {'mumap','ct','mu-map'}
                    simSet.CTscanNum  = i;
                case {'pim','parametric image','parametricimage','parametric_image'}
                    simSet.PIMscanNum = i;
                case {'pt','pet','pet4d','pt pristine','pristine','image'}
                    simSet.PTscanNum  = i;
                otherwise
                    %% Reset pointer
                    handles.buttonRun.Value = 0;
                    handles.running.String  = '';
                    set(handles.figure1, 'pointer', oldpointer)
                    
                    %% Break
                    fprintf('Unrecognized data type "%s". Aborting.\n',types{i});
                    return;
            end
        end
    end
    
    %% Check that data input is OK
    if handles.inputPIM.Value == 1
        if ~isfield(simSet,'PIMscanNum') 
            %% Break
            fprintf('Missing parametric image data in structure "%s". Aborting.\n',handles.dataInputName.String);
            return
        end
    end
    if handles.inputDYN.Value == 1
        if ~isfield(simSet,'PTscanNum')            
            %% Break
            fprintf('Missing dynamic image data in structure "%s". Aborting.\n',handles.dataInputName.String);
            return
        end
    end
    
    %% Add PT scan num if needed
    if handles.inputPIM.Value ==1
        simSet.PTscanNum = numel(types)+1;
    end
    
    %% Count scale factor
    scaleFactor = handles.scaleFactor.Value;
    
    %% Run simulation
    try
        clock = tic;
        if handles.inputPIM.Value==1 %Parametric image input
            [data,FBP,OSEM,OSEMpsf,counts,countsNoise,nFWprompts,FWtrues,FWscatters,FWrandoms,wcc] = ...
                gui_Dynamic_main_kineticModelling(data,simSet,scaleFactor);
        else %Dynamic image input
            [FBP,OSEM,OSEMpsf,counts,countsNoise,nFWprompts,FWtrues,FWscatters,FWrandoms,wcc] = ...
                gui_Dynamic_main_dynamicImage(data,simSet,scaleFactor);
        end
        timing = toc(clock);
    catch ME
        rethrow(ME)
    end
    
    %% Save to workspace
    fprintf('Transferring results to workspace... ')
    editBoxes = findobj(allchild(handles.uipanel6),'style','edit');
    for i=1:numel(editBoxes)
        if strcmp(handles.(editBoxes(i).Tag).Enable,'on')
           assignin('base',editBoxes(i).String,eval(editBoxes(i).Tag(1:end-4)))
        end
    end
    fprintf('Done!\n')
    
    %% Indicate stop running with pointer
    handles.buttonRun.Value = 0;
    set(handles.figure1,'pointer',oldpointer);
    handles.running.String  = '';
    
    fprintf('\nAll done with simulation (%.0f sec | %.2f min).\n',timing,timing/60)
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [hObject,error] = checkIfNumeric(hObject)
    error = 0;
    if ( isempty(hObject.Value) && ~isempty(hObject.String) ) %non-numeric, non empty entry
        fprintf('Value "%s" in box "%s" is not a number.\n',get(hObject,'String'),get(hObject,'Tag'));
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    else
        set(hObject,'BackgroundColor',[1 1 1]);     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hObject,error] = checkIfAboveZero(hObject)
    error = 0;
    if hObject.Value<=0
        fprintf('Value "%s" in box "%s" has to be larger than zero.\n',get(hObject,'Value'),get(hObject,'Tag'));
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    else
        set(hObject,'BackgroundColor',[1 1 1]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hObject,error] = checkIfNonNeg(hObject)
    error = 0;
    if hObject.Value<0
        fprintf('Value "%s" in box "%s" has to be a positive number.\n',get(hObject,'Value'),get(hObject,'Tag'));
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    else
        set(hObject,'BackgroundColor',[1 1 1]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCleanupFun(handles,oldpointer)
    diary off
    handles.buttonRun.Value  = 0;
    handles.running.String   = '';
    set(handles.figure1,'pointer',oldpointer);
    fprintf('\nAborting run!\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
