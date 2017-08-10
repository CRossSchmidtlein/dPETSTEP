function varargout = dPETSTEPgui_fit(varargin)
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

    % Last Modified by GUIDE v2.5 08-Aug-2017 16:57:00

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @dPETSTEPgui_fit_OpeningFcn, ...
        'gui_OutputFcn',  @dPETSTEPgui_fit_OutputFcn, ...
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
function dPETSTEPgui_fit_OpeningFcn(hObject, eventdata, handles, varargin)
    
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
            panelGroups = {'FitType','dynamicPETName'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
        elseif strcmp(panel{i},'uipanel3')
            panelGroups = {'inputFunction','inputFunctionROI'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
            uistack(handles.inputFuncExt,'bottom')
            uistack(handles.inputFuncID,'bottom')
        elseif strcmp(panel{i},'uipanel4')
            panelGroups = {'StopCriteria','Solver','Interp','frameWeight'};
            for j = 1:numel(panelGroups)
                eval(['uistack( findobj(handles.'  panel{i} ...
                    ', ''Tag'',''' panelGroups{j} '''),''bottom'')'])
            end
        end
    end

    %% Update handles structure
    guidata(hObject, handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = dPETSTEPgui_fit_OutputFcn(hObject, eventdata, handles) 

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataFromWS_Callback(hObject, eventdata, handles)
    
    %% Reset color
    set(handles.dataName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function dataFromFile_Callback(hObject, eventdata, handles)

    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat'},'Pick a file');
    
    %% Reset color
    set(handles.dataName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.dataName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function dataName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function voxelwiseFit_Callback(hObject, eventdata, handles)
    %% Disable ROI info
    handles.roiNameText.Enable = 'off';
    handles.roiName.Enable     = 'off';
    handles.roiFromWS.Enable   = 'off';
    handles.roiFromFile.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function roiFit_Callback(hObject, eventdata, handles)

    %% Disable ROI info
    handles.roiNameText.Enable = 'on';
    handles.roiName.Enable     = 'on';
    handles.roiFromWS.Enable   = 'on';
    handles.roiFromFile.Enable = 'on';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function roiFromWS_Callback(hObject, eventdata, handles)
    %% Reset color
    set(handles.roiName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function roiFromFile_Callback(hObject, eventdata, handles)
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat'},'Pick a file');
    
    %% Reset color
    set(handles.roiName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.roiName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function roiName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OneTissue_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable     = 'off';
    handles.textNumExp.Enable = 'off';
      
    %% Enable relevant parameter settings
    handles.textVa.Enable    = 'on';
    handles.includeVa.Enable = 'on';
        includeVa_Callback(handles.includeVa,[],handles);
    handles.k1Init.Enable    = 'on';
    handles.k1Lower.Enable   = 'on';
    handles.k1Upper.Enable   = 'on';
    handles.k2Init.Enable    = 'on';
    handles.k2Lower.Enable   = 'on';
    handles.k2Upper.Enable   = 'on';
    handles.textK1.Enable    = 'on';
    handles.textk2.Enable    = 'on';
    
    %% Disable irrelevant parameter settings
    handles.textk3.Enable    = 'off';
    handles.textk4.Enable    = 'off';
    handles.textR1.Enable    = 'off';
    handles.textBPnd.Enable  = 'off';
    handles.textP.Enable     = 'off';
    handles.textPinit.Enable  = 'off';
    handles.textPlower.Enable = 'off';
    handles.textPupper.Enable = 'off';
    handles.textPuseSep.Enable = 'off';
    handles.k3Init.Enable    = 'off';
    handles.k3Lower.Enable   = 'off';
    handles.k3Upper.Enable   = 'off';
    handles.k4Init.Enable    = 'off';
    handles.k4Lower.Enable   = 'off';
    handles.k4Upper.Enable   = 'off';
    handles.r1Init.Enable    = 'off';
    handles.r1Lower.Enable   = 'off';
    handles.r1Upper.Enable   = 'off';
    handles.bpndInit.Enable  = 'off';
    handles.bpndLower.Enable = 'off';
    handles.bpndUpper.Enable = 'off';
    handles.pInit.Enable     = 'off';
    handles.pLower.Enable    = 'off';
    handles.pUpper.Enable    = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function TwoTissue_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable     = 'off';
    handles.textNumExp.Enable = 'off';
       
    %% Enable relevant parameter settings
    handles.textVa.Enable    = 'on';
    handles.includeVa.Enable = 'on';
        includeVa_Callback(handles.includeVa,[],handles);
    handles.textK1.Enable    = 'on';
    handles.textk2.Enable    = 'on';
    handles.textk3.Enable    = 'on';
    handles.textk4.Enable    = 'on';
    handles.k1Init.Enable    = 'on';
    handles.k1Lower.Enable   = 'on';
    handles.k1Upper.Enable   = 'on';
    handles.k2Init.Enable    = 'on';
    handles.k2Lower.Enable   = 'on';
    handles.k2Upper.Enable   = 'on';
    handles.k3Init.Enable    = 'on';
    handles.k3Lower.Enable   = 'on';
    handles.k3Upper.Enable   = 'on';
    handles.k4Init.Enable    = 'on';
    handles.k4Lower.Enable   = 'on';
    handles.k4Upper.Enable   = 'on';
    
    %% Disable irrelevant parameter settings
    handles.textR1.Enable    = 'off';
    handles.textBPnd.Enable  = 'off';
    handles.textP.Enable     = 'off';
    handles.textPinit.Enable  = 'off';
    handles.textPlower.Enable = 'off';
    handles.textPupper.Enable = 'off';
    handles.textPuseSep.Enable = 'off';
    handles.r1Init.Enable    = 'off';
    handles.r1Lower.Enable   = 'off';
    handles.r1Upper.Enable   = 'off';
    handles.bpndInit.Enable  = 'off';
    handles.bpndLower.Enable = 'off';
    handles.bpndUpper.Enable = 'off';
    handles.pInit.Enable     = 'off';
    handles.pLower.Enable    = 'off';
    handles.pUpper.Enable    = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function FRTM_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable     = 'off';
    handles.textNumExp.Enable = 'off';
    
    %% Enable relevant parameter settings
    handles.r1Init.Enable    = 'on';
    handles.r1Lower.Enable   = 'on';
    handles.r1Upper.Enable   = 'on';
    handles.bpndInit.Enable  = 'on';
    handles.bpndLower.Enable = 'on';
    handles.bpndUpper.Enable = 'on';
    handles.textR1.Enable    = 'on';
    handles.textBPnd.Enable  = 'on';
    
    %% Disable irrelevant parameter settings
    handles.includeVa.Enable = 'off';
    handles.textVa.Enable    = 'off';
    handles.textP.Enable      = 'off';
    handles.textPinit.Enable  = 'off';
    handles.textPlower.Enable = 'off';
    handles.textPupper.Enable = 'off';
    handles.textPuseSep.Enable = 'off';
    handles.vaInit.Enable    = 'off';
    handles.vaLower.Enable   = 'off';
    handles.vaUpper.Enable   = 'off';
    handles.textK1.Enable    = 'off';
    handles.textk2.Enable    = 'off';
    handles.textk3.Enable    = 'off';
    handles.textk4.Enable    = 'off';
    handles.textP.Enable     = 'off';
    handles.k1Init.Enable    = 'off';
    handles.k1Lower.Enable   = 'off';
    handles.k1Upper.Enable   = 'off';
    handles.k2Init.Enable    = 'off';
    handles.k2Lower.Enable   = 'off';
    handles.k2Upper.Enable   = 'off';
    handles.k3Init.Enable    = 'off';
    handles.k3Lower.Enable   = 'off';
    handles.k3Upper.Enable   = 'off';
    handles.k4Init.Enable    = 'off';
    handles.k4Lower.Enable   = 'off';
    handles.k4Upper.Enable   = 'off';
    handles.pInit.Enable     = 'off';
    handles.pLower.Enable    = 'off';
    handles.pUpper.Enable    = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function SRTM_Callback(hObject, eventdata, handles)

    %% Disable number of exponentials box
    handles.numExp.Enable     = 'off';
    handles.textNumExp.Enable = 'off';
      
    %% Enable relevant parameter settings
    handles.textR1.Enable    = 'on';
    handles.textBPnd.Enable  = 'on';
    handles.r1Init.Enable    = 'on';
    handles.r1Lower.Enable   = 'on';
    handles.r1Upper.Enable   = 'on';
    handles.bpndInit.Enable  = 'on';
    handles.bpndLower.Enable = 'on';
    handles.bpndUpper.Enable = 'on';
        
    %% Disable irrelevant parameter settings
    handles.textK1.Enable    = 'off';
    handles.textk2.Enable    = 'off';
    handles.textk3.Enable    = 'off';
    handles.textk4.Enable    = 'off';
    handles.textP.Enable      = 'off';
    handles.textPinit.Enable  = 'off';
    handles.textPlower.Enable = 'off';
    handles.textPupper.Enable = 'off';
    handles.textPuseSep.Enable = 'off';
    handles.k1Init.Enable    = 'off';
    handles.k1Lower.Enable   = 'off';
    handles.k1Upper.Enable   = 'off';
    handles.k2Init.Enable    = 'off';
    handles.k2Lower.Enable   = 'off';
    handles.k2Upper.Enable   = 'off';
    handles.k3Init.Enable    = 'off';
    handles.k3Lower.Enable   = 'off';
    handles.k3Upper.Enable   = 'off';
    handles.k4Init.Enable    = 'off';
    handles.k4Lower.Enable   = 'off';
    handles.k4Upper.Enable   = 'off';
    handles.includeVa.Enable = 'off';
    handles.textVa.Enable    = 'off';
    handles.vaInit.Enable    = 'off';
    handles.vaLower.Enable   = 'off';
    handles.vaUpper.Enable   = 'off';
    handles.pInit.Enable     = 'off';
    handles.pLower.Enable    = 'off';
    handles.pUpper.Enable    = 'off';  
    
    %% Update handles structure
    guidata(hObject, handles);
end

function SumExp_Callback(hObject, eventdata, handles)

    %% Enable number of exponentials box
    handles.numExp.Enable     = 'on';
    handles.textNumExp.Enable = 'on';
    
    %% Enable relevant parameter settings
    handles.vaUpper.Enable   = 'on';
    handles.textP.Enable     = 'on';
    handles.textPinit.Enable  = 'on';
    handles.textPlower.Enable = 'on';
    handles.textPupper.Enable = 'on';
    handles.textPuseSep.Enable = 'on';
    handles.pInit.Enable     = 'on';
    handles.pLower.Enable    = 'on';
    handles.pUpper.Enable    = 'on';
    handles.includeVa.Enable = 'on';
    handles.textVa.Enable    = 'on';
        includeVa_Callback(handles.includeVa,[],handles);
     
    %% Disable irrelevant parameter settings
    handles.textK1.Enable    = 'off';
    handles.textk2.Enable    = 'off';
    handles.textk3.Enable    = 'off';
    handles.textk4.Enable    = 'off';
    handles.textR1.Enable    = 'off';
    handles.textBPnd.Enable  = 'off';
    handles.k1Init.Enable    = 'off';
    handles.k1Lower.Enable   = 'off';
    handles.k1Upper.Enable   = 'off';
    handles.k2Init.Enable    = 'off';
    handles.k2Lower.Enable   = 'off';
    handles.k2Upper.Enable   = 'off';
    handles.k3Init.Enable    = 'off';
    handles.k3Lower.Enable   = 'off';
    handles.k3Upper.Enable   = 'off';
    handles.k4Init.Enable    = 'off';
    handles.k4Lower.Enable   = 'off';
    handles.k4Upper.Enable   = 'off';
    handles.r1Init.Enable    = 'off';
    handles.r1Lower.Enable   = 'off';
    handles.r1Upper.Enable   = 'off';
    handles.bpndInit.Enable  = 'off';
    handles.bpndLower.Enable = 'off';
    handles.bpndUpper.Enable = 'off';
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function frameName_Callback(hObject, eventdata, handles)

    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inputFuncID_Callback(hObject, eventdata, handles)

    %% Disable external info
    handles.inputFuncNameText.Enable = 'off';
    handles.inputFuncNameUnit.Enable = 'off';
    handles.inputFuncName.Enable     = 'off';
    handles.ifFromWS.Enable          = 'off';
    handles.ifFromFile.Enable        = 'off';
    
    %% Enable ROI info
    handles.inputFuncROINameText.Enable = 'on';
    handles.inputFuncROIName.Enable     = 'on';
    handles.ifROIFromWS.Enable          = 'on';
    handles.ifROIFromFile.Enable        = 'on';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function inputFuncExt_Callback(hObject, eventdata, handles)
    %% Disable external info
    handles.inputFuncNameText.Enable = 'on';
    handles.inputFuncNameUnit.Enable = 'on';
    handles.inputFuncName.Enable     = 'on';
    handles.ifFromWS.Enable          = 'on';
    handles.ifFromFile.Enable        = 'on';
    
    %% Enable ROI info
    handles.inputFuncROINameText.Enable = 'off';
    handles.inputFuncROIName.Enable     = 'off';
    handles.ifROIFromWS.Enable          = 'off';
    handles.ifROIFromFile.Enable        = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function ifROIFromWS_Callback(hObject, eventdata, handles)
    %% Reset color
    set(handles.inputFuncName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function ifROIFromFile_Callback(hObject, eventdata, handles)
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat'},'Pick a file');
    
    %% Reset color
    set(handles.inputFuncName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.inputFuncName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function inputFuncROIName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function ifFromWS_Callback(hObject, eventdata, handles)
    %% Update handles structure
    guidata(hObject, handles);
end

function ifFromFile_Callback(hObject, eventdata, handles)
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat;*.txt;*.csv*'},'Pick a file');
    
    %% Put file name in box
    handles.inputFuncName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
    
end

function inputFuncName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ones_Callback(hObject, eventdata, handles)
    
    %% Disable halflife
    handles.halflife.Enable     = 'off';
    handles.halflifeText.Enable = 'off';
    handles.halflifeUnit.Enable = 'off';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function w1_Callback(hObject, eventdata, handles)

    %% Enable halflife
    handles.halflife.Enable     = 'on';
    handles.halflifeText.Enable = 'on';
    handles.halflifeUnit.Enable = 'on';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function w2_Callback(hObject, eventdata, handles)
    %% Enable halflife
    handles.halflife.Enable     = 'on';
    handles.halflifeText.Enable = 'on';
    handles.halflifeUnit.Enable = 'on';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function w3_Callback(hObject, eventdata, handles)
    %% Enable halflife
    handles.halflife.Enable     = 'on';
    handles.halflifeText.Enable = 'on';
    handles.halflifeUnit.Enable = 'on';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function w4_Callback(hObject, eventdata, handles)
    %% Enable halflife
    handles.halflife.Enable     = 'on';
    handles.halflifeText.Enable = 'on';
    handles.halflifeUnit.Enable = 'on';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function w6_Callback(hObject, eventdata, handles)
    %% Enable halflife
    handles.halflife.Enable     = 'on';
    handles.halflifeText.Enable = 'on';
    handles.halflifeUnit.Enable = 'on';
    
    %% Disable read data
    handles.weightName.Enable     = 'off';
    handles.weightNameText.Enable = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

function weightFromWS_Callback(hObject, eventdata, handles)
    %% Enable read data
    handles.weightName.Enable     = 'on';
    handles.weightNameText.Enable = 'on';
    
    %% Reset color
    set(handles.weightName,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function weightFromFile_Callback(hObject, eventdata, handles)
    %% Enable read data
    handles.weightName.Enable     = 'on';
    handles.weightNameText.Enable = 'on';
    
    %% Pick file
    [fileName,pathName] = uigetfile({'*.mat;*.txt;*.csv*'},'Pick a file');
    
    %% Reset color
    set(handles.weightName,'BackgroundColor', [1 1 1] );
    
    %% Put file name in box
    handles.weightName.String = fullfile(pathName,fileName);
    
    %% Update handles structure
    guidata(hObject, handles);
end

function weightName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );

    %% Update handles structure
    guidata(hObject, handles);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trustregionreflective_Callback(hObject, eventdata, handles)
    %% Enable bounds
    model = get(get(handles.KineticModel,'SelectedObject'), 'String');
    if ismember(model,{'1-Tissue','2-Tissue'})
        handles.k1Lower.Enable   = 'on';
        handles.k1Upper.Enable   = 'on';
        handles.k2Lower.Enable   = 'on';
        handles.k2Upper.Enable   = 'on';
    end
    if ismember(model,{'2-Tissue'})
        handles.k3Lower.Enable   = 'on';
        handles.k3Upper.Enable   = 'on';
        handles.k4Lower.Enable   = 'on';
        handles.k4Upper.Enable   = 'on';
    end
    if ismember(model,{'1-Tissue','2-Tissue','SumExp'})
        if handles.includeVa.Value == 1
            handles.vaLower.Enable   = 'on';
            handles.vaUpper.Enable   = 'on';
        end
    end
    if ismember(model,{'FRTM','SRTM'})
        handles.r1Lower.Enable   = 'on';
        handles.r1Upper.Enable   = 'on';
        handles.bpndLower.Enable = 'on';
        handles.bpndUpper.Enable = 'on';
    end
    if ismember(model,{'SumExp'})
        handles.pLower.Enable    = 'on';
        handles.pUpper.Enable    = 'on';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function levenbergmarquardt_Callback(hObject, eventdata, handles)

    %% Disable bounds
    handles.k1Lower.Enable   = 'off';
    handles.k1Upper.Enable   = 'off';
    handles.k2Lower.Enable   = 'off';
    handles.k2Upper.Enable   = 'off';
    handles.k3Lower.Enable   = 'off';
    handles.k3Upper.Enable   = 'off';
    handles.k4Lower.Enable   = 'off';
    handles.k4Upper.Enable   = 'off';
    handles.vaLower.Enable   = 'off';
    handles.vaUpper.Enable   = 'off';
    handles.r1Lower.Enable   = 'off';
    handles.r1Upper.Enable   = 'off';
    handles.bpndLower.Enable = 'off';
    handles.bpndUpper.Enable = 'off';
    handles.pLower.Enable    = 'off';
    handles.pUpper.Enable    = 'off';
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MaxIterCheck_Callback(hObject, eventdata, handles)

    %% Disable MaxIter box or not
    if hObject.Value==1 %Box ticked
        handles.MaxIter.Enable = 'on';
        handles.MaxIter.String = '';
    else %Box unticked
        handles.MaxIter.Enable = 'off';
        handles.MaxIter.String = '';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function MaxIter_Callback(hObject, eventdata, handles)
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

function TolFunCheck_Callback(hObject, eventdata, handles)
    %% Disable MaxIter box or not
    if hObject.Value==1 %Box ticked
        handles.TolFun.Enable = 'on';
        handles.TolFun.String = '';
    else %Box unticked
        handles.TolFun.Enable = 'off';
        handles.TolFun.String = '';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function TolFun_Callback(hObject, eventdata, handles)
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

function TolXCheck_Callback(hObject, eventdata, handles)
    %% Disable MaxIter box or not
    if hObject.Value==1 %Box ticked
        handles.TolX.Enable = 'on';
        handles.TolX.String = '';
    else %Box unticked
        handles.TolX.Enable = 'off';
        handles.TolX.String = '';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end
function TolX_Callback(hObject, eventdata, handles)
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

function noCPU_Callback(hObject, eventdata, handles)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unitSec_Callback(hObject, eventdata, handles)
    %% Update handles structure
    guidata(hObject, handles);
end

function unitMin_Callback(hObject, eventdata, handles)
    %% Update handles structure
    guidata(hObject, handles);
end

function k1Init_Callback(hObject, eventdata, handles)
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

function k2Init_Callback(hObject, eventdata, handles)
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

function k3Init_Callback(hObject, eventdata, handles)
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

function k4Init_Callback(hObject, eventdata, handles)
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

function vaInit_Callback(hObject, eventdata, handles)
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

function r1Init_Callback(hObject, eventdata, handles)
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

function bpndInit_Callback(hObject, eventdata, handles)
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

function k1Upper_Callback(hObject, eventdata, handles)
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

function k2Upper_Callback(hObject, eventdata, handles)
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

function k3Upper_Callback(hObject, eventdata, handles)
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

function k4Upper_Callback(hObject, eventdata, handles)
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

function r1Upper_Callback(hObject, eventdata, handles)
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

function bpndUpper_Callback(hObject, eventdata, handles)
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

function vaUpper_Callback(hObject, eventdata, handles)
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

function k1Lower_Callback(hObject, eventdata, handles)
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

function k2Lower_Callback(hObject, eventdata, handles)
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

function k3Lower_Callback(hObject, eventdata, handles)
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

function k4Lower_Callback(hObject, eventdata, handles)
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

function r1Lower_Callback(hObject, eventdata, handles)
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

function bpndLower_Callback(hObject, eventdata, handles)
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

function vaLower_Callback(hObject, eventdata, handles)
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

function includeVa_Callback(hObject, eventdata, handles)
    %% Enable parameter settings
    if hObject.Value==1
%         handles.textVa.Enable  = 'on';
        handles.vaInit.Enable  = 'on';
        handles.vaUpper.Enable = 'on';
        handles.vaLower.Enable = 'on';
    else
%         handles.textVa.Enable  = 'off';
        handles.vaInit.Enable  = 'off';
        handles.vaUpper.Enable = 'off';
        handles.vaLower.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end


function pInit_Callback(hObject, eventdata, handles)   

    if ~isempty(eventdata)
        %% Check if numeric value
        hObject = checkIfNumericCommaSep(hObject);
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

function pLower_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Check if numeric value
        hObject = checkIfNumericCommaSep(hObject);
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

function pUpper_Callback(hObject, eventdata, handles)

    if ~isempty(eventdata)
        %% Check if numeric value
        hObject = checkIfNumericCommaSep(hObject);
        
        %% Update handles structure
        guidata(hObject, handles);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pim_Callback(hObject, eventdata, handles)
    %% Enable output
    if hObject.Value==1
        handles.pimName.Enable = 'on';
    else
        handles.pimName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function pimName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end

function fitInfo_Callback(hObject, eventdata, handles)
    %% Enable output
    if hObject.Value==1
        handles.fitInfoName.Enable = 'on';
    else
        handles.fitInfoName.Enable = 'off';
    end
    
    %% Update handles structure
    guidata(hObject, handles);
end

function fitInfoName_Callback(hObject, eventdata, handles)
    %% Background to white
    set(hObject,'BackgroundColor', [1 1 1] );
    
    %% Update handles structure
    guidata(hObject, handles);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
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
    handles = gui_loadSettings(handles,'fit');
    
    %% Update handles structure
    guidata(hObject, handles); 
end


function buttonReset_Callback(hObject, eventdata, handles)
    
    %% Reset all values
    % Reset all backgrounds to white
    editBoxes = findobj(gcf,'style','edit');
    for i=1:numel(editBoxes)
        handles.(editBoxes(i).Tag).BackgroundColor = [1 1 1];
    end
    
    % Radio buttons
    set(handles.dataFromWS,'Value',1);
    set(handles.frameFromWS,'Value',1);
    set(handles.voxelwiseFit,'Value',1);
        voxelwiseFit_Callback(handles.voxelwiseFit,[], handles)
    set(handles.inputFuncID,'Value',1);
        inputFuncID_Callback(handles.inputFuncID,[], handles)
    set(handles.OneTissue,'Value',1);
        OneTissue_Callback(handles.OneTissue,[], handles)
    set(handles.ones,'Value',1);
        ones_Callback(handles.ones,[], handles)
    set(handles.linear,'Value',1);
        linear_Callback(handles.linear,[], handles)
    set(handles.trustregionreflective,'Value',1);
    set(handles.unitSec,'Value',1);
    
    % Edit boxes
    editBoxes = findobj(gcf,'style','edit');
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

    % Checkboxes
    checkBoxes = findobj(allchild(handles.uipanel5),'style','checkbox');
    for i=1:numel(checkBoxes)
        set(handles.(checkBoxes(i).Tag),'Value',1);
    end
    set(handles.MaxIterCheck,'Value',0);
        MaxIterCheck_Callback(handles.MaxIterCheck,[], handles)
    set(handles.TolFunCheck,'Value',0);
        TolFunCheck_Callback(handles.TolFunCheck,[], handles)
    set(handles.TolXCheck,'Value',0);
        TolXCheck_Callback(handles.TolXCheck,[], handles)
    set(handles.pim,'Value',1);
        pim_Callback(handles.pim,[], handles)
    set(handles.fitInfo,'Value',0);
        fitInfo_Callback(handles.fitInfo,[], handles)
        
    % Special cases
    set(handles.timeStep,'Value',1);          set(handles.timeStep,'String','1');
    set(handles.noCPU,'Value',1);             set(handles.noCPU,'String','1');
    set(handles.dataName,'Value',[]);         set(handles.dataName,'String','image4D');
    set(handles.roiName,'Value',[]);          set(handles.roiName,'String','ROI');
    set(handles.frameName,'Value',[]);        set(handles.frameName,'String','midFrame');
    set(handles.weightName,'Value',[]);       set(handles.weightName,'String','weight');
    set(handles.inputFuncROIName,'Value',[]); set(handles.inputFuncROIName,'String','ROI_Cif');
    set(handles.inputFuncName,'Value',[]);    set(handles.inputFuncName,'String','Cif');
    set(handles.pimName,'Value',[]);          set(handles.pimName,'String','pFitted');
    set(handles.fitInfoName,'Value',[]);      set(handles.fitInfoName,'String','fitInfo');
   
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
    varName    = {handles.dataName   handles.roiName   handles.frameName   handles.inputFuncROIName handles.inputFuncName handles.weightName};
    toggleName = {handles.dataFromWS handles.roiFromWS handles.frameFromWS handles.ifROIFromWS      handles.ifFromWS      handles.weightFromWS}; 
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
                if strcmp(ext,'.mat') % MAT-file
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
    boxes = findobj(gcf,'style','edit');
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
       
    %% Check that correct # of parameters is sum of exponentials
    if handles.SumExp.Value == 1 %chosen
         noExp = handles.numExp.Value;
         noP    = noExp*2; %number of parameters
         name   = {'pInit','pLower','pUpper'};
         for i=1:numel(name)
             num = numel( strsplit( get(handles.(name{i}),'String'),',' ) );
             if num ~= noP
                 fprintf('Wrong number of parameters in box "%s". Should be %d.\n',get(handles.(name{i}),'Tag'),noP);
                 set(handles.(name{i}),'BackgroundColor',[1 0.4 0.4]);
                 error = 1;
             end
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
 
    %% Create fitting settings
    options           = optimset('Display','off');
    options.Algorithm = get(get(handles.Solver,'SelectedObject'), 'String');
    boxes             = findobj(handles.StopCriteria,'style','edit');
    for i=1:numel(boxes)
        if strcmp(boxes(i).Enable,'on')
            options.(boxes(i).Tag) = boxes(i).Value;
        end
    end
    
    %% Read 4D image from WS or file
    if handles.dataFromWS.Value==1 %selected read from WS
        image4D   = evalin('base',handles.dataName.String);
    else
        [~,~,ext] = fileparts(handles.dataName.String);
        if strcmp(ext,'.mat') %MAT-file
            image4D = cell2mat(struct2cell( load(handles.dataName.String) ));
        else %Text file
            image4D = load(handles.dataName.String);
        end
    end

    %% Read 4D image from WS or file
    if handles.frameFromWS.Value==1 %selected read from WS
        midFrame   = evalin('base',handles.frameName.String);
    else
        [~,~,ext] = fileparts(handles.frameName.String);
        if strcmp(ext,'.mat') %MAT-file
            midFrame = cell2mat(struct2cell( load(handles.frameName.String) ));
        else %Text file
            midFrame = load(handles.frameName.String);
        end
    end
      
    %% Input function
    [Cif,CifMask] = getInputFunction(handles);
    
     %% Make sure input function and dynamic data have equally many time points
    nt = size(image4D,4);
    if numel(Cif)~=nt 
        fprintf('Time dimension mismatch! Input function has %d time points, dynamic data has %d.\n',numel(Cif),nt);
        set(handles.inputFuncName,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    end
    if numel(midFrame)~=nt 
        fprintf('Time dimension mismatch! Mid frame has %d time points, dynamic data has %d.\n',numel(midFrame),nt);
        set(handles.frameName,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    end
    % Break if error
    if error; return; end
    
    %% Frame weight
    if handles.weightFromWS.Value==1 %selected read from WS
        w   = evalin('base',handles.weightName.String);
    elseif handles.weightFromFile.Value==1
        [~,~,ext] = fileparts(handles.weightName.String);
        if strcmp(ext,'.mat') %MAT-file
            w = cell2mat(struct2cell( load(handles.weightName.String) ));
        else %Text file
            w = load(handles.weightName.String);
        end
    else
        w   = get(get(handles.frameWeight,'SelectedObject'), 'Tag');
    end
    
    %% Initial value and lower+upper bounds
    p0 = []; lb = []; ub = [];
    order    = {'k1','k2','k3','k4','r1','bpnd','p','va'};
    boxes    = findobj(handles.uipanel5,'style','edit');
    for i=1:numel(boxes)
        if strcmp(boxes(i).Enable,'on')
            if strcmp(boxes(i).Tag(1),'p')
                % Cell to number
                numCell = strsplit( get(boxes(i),'String'),',' );
                num     = zeros(numel(numCell),1);
                for j=1:numel(numCell)
                    num(j) = str2num(numCell{j});
                end
                % Initial value, lower or upper bound
                if strcmp(boxes(i).Tag(end-3:end),'Init')
                    p0 = [num (1:numel(num))'];
                elseif strcmp(boxes(i).Tag(end-4:end),'Lower')
                    lb = [num (1:numel(num))'];
                elseif strcmp(boxes(i).Tag(end-4:end),'Upper')
                    ub = [num (1:numel(num))'];
                end
                        
            elseif strcmp(boxes(i).Tag(end-3:end),'Init')
                param = boxes(i).Tag(1:end-4);
                [~,n] = ismember(param,order);
                p0 = [p0 ; boxes(i).Value n];
            elseif strcmp(boxes(i).Tag(end-4:end),'Lower')
                param = boxes(i).Tag(1:end-5);
                [~,n] = ismember(param,order);
                lb = [lb ; boxes(i).Value n];
            elseif strcmp(boxes(i).Tag(end-4:end),'Upper')
                param = boxes(i).Tag(1:end-5);
                [~,n] = ismember(param,order);
                ub = [ub ; boxes(i).Value n];
            end
        end
    end
    % Use appropriate unit (default 1/sec). 
    if handles.unitMin.Value == 1 %Convert 1/min to 1/sec
        % Only for k1, k2, k3, k4 and p.
        rateConst = [1 2 3 4 7]; %indices from cell "order" above. 
        for i=1:numel(rateConst)
            [yes,ind] = ismember(rateConst(i),p0(:,2));
            if yes
                p0(ind,1) = p0(ind,1)/2;
            end
            [yes,ind] = ismember(rateConst(i),lb(:,2));
            if yes
                lb(ind,1) = lb(ind,1)/2;
            end
            [yes,ind] = ismember(rateConst(i),ub(:,2));
            if yes
                ub(ind,1) = ub(ind,1)/2;
            end
        end
    end
    % Sort in right order
    [~,ind] = sort(p0(:,2),'ascend');     P0 = p0(ind,1);
    if ~isempty(lb)
        [~,ind] = sort(lb(:,2),'ascend'); LB = lb(ind,1);
    end
    if ~isempty(ub)
        [~,ind] = sort(ub(:,2),'ascend'); UB = ub(ind,1);
    end
    
    %% Chosen kinetic model
    model = get(get(handles.KineticModel,'SelectedObject'), 'String');
    
    %% ROI fitting (maybe)
    ROIMask = [];
    if handles.roiFit.Value==1
        if handles.roiFromWS.Value==1 %selected read from WS
            ROIMask  = evalin('base',handles.roiName.String);
        else
            [~,~,ext] = fileparts(handles.roiName.String);
            if strcmp(ext,'.mat') %MAT-file
                ROIMask = cell2mat(struct2cell( load(handles.roiName.String) ));
            else %Text file
                ROIMask = load(handles.roiName.String);
            end
        end
    end

    %% Halflife
    halflife = [];
    if strcmp(handles.halflife.Enable,'on')
        halflife = handles.halflife.Value;
    end
      
    %% Number of CPUs
    noCPU = handles.noCPU.Value;
    
    %% Interpolation and dt
    interpMethod = get(get(handles.Interp,'SelectedObject'), 'String');
    dt           = handles.timeStep.Value;
    
    %% Input string of arguments
    inputStr = ['''image'',image4D,''w'',w,''model'',model,''midFrame'',midFrame,'...
        '''p0'',P0,''interpMethod'',interpMethod,''dt'',dt,''options'',options,''noCPU'',noCPU'];
    
    %% Optional arguments
    if ~isempty(Cif)
        if ismember(model,{'1-Tissue','2-Tissue','SumExp'});
            inputStr = [inputStr ',''Cp'',Cif'];
        else
            inputStr = [inputStr ',''Cref'',Cif'];
        end
    end
    if ~isempty(CifMask)
        if ismember(model,{'1-Tissue','2-Tissue','SumExp'});
            inputStr = [inputStr ',''CpMask'',CifMask'];
        else
            inputStr = [inputStr ',''CrefMask'',CifMask'];
        end
    end
    if ~isempty(ROIMask)
        inputStr = [inputStr ',''ROIMask'',ROIMask'];
    end
    if ~isempty(halflife)
        inputStr = [inputStr ',''halflife'',halflife'];
    end
    if handles.SumExp.Value==1
        noExp   = handles.numExp.Value;
        inputStr = [inputStr ',''noExp'',noExp'];
    end
    if handles.trustregionreflective.Value==1
        inputStr = [inputStr ',''lowerBound'',LB,''upperBound'',UB'];
    end
    
    %% Run fitting 
    try
        clock  = tic;
        eval( ['[pimFitted,fitInfo] = modelFitting_main(' inputStr ');'] );
        timing = toc(clock);
    catch ME
        rethrow(ME)
    end
    
    %% Save to workspace
    if handles.pim.Value==1 || handles.fitInfo.Value==1; 
        fprintf('Transferring results to workspace... ')
    end
    if handles.pim.Value==1
        assignin('base',handles.pimName.String,pimFitted)
    end
    if handles.fitInfo.Value==1
        assignin('base',handles.fitInfoName.String,fitInfo)
    end
    if handles.pim.Value==1 || handles.fitInfo.Value==1;
        fprintf('Done!\n')
    end
    
    %% Indicate stop running with pointer and text
    handles.buttonRun.Value = 0;
    set(handles.figure1,'pointer',oldpointer);
    handles.running.String  = '';
    
    fprintf('\nAll done with fitting (%.0f sec | %.2f min).\n',timing,timing/60)
    
    %% Update handles structure
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hObject,error] = checkIfNumeric(hObject)
    error = 0;
    if isempty(hObject.Value) && ~isempty(hObject.String) %non-numeric
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
        fprintf('Value "%s" in box "%s" has to be positive number.\n',get(hObject,'Value'),get(hObject,'Tag'));
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    else
        set(hObject,'BackgroundColor',[1 1 1]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hObject,error] = checkIfNumericCommaSep(hObject)
    error = 0;
    num   = strsplit( get(hObject,'String'),',' )';

    %% Check if numeric value
     for i=1:numel(num)
        [n,ok] = str2num(num{i});
        if ok==1 && (~isreal(n) || n<0)
            ok = 0;
        end
    end
    if ~ok %non-numeric
        fprintf('Values "%s" in box "%s" are not comma separated, positive numbers.\n',get(hObject,'String'),get(hObject,'Tag'));
        set(hObject,'BackgroundColor',[1 0.4 0.4]);
        error = 1;
    else
        set(hObject,'BackgroundColor',[1 1 1]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCleanupFun(handles,oldpointer)
    handles.buttonRun.Value  = 0;
    handles.running.String   = '';
    set(handles.figure1,'pointer',oldpointer);
    fprintf('\nAborting run!\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%