function varargout = dPETSTEP(varargin)
    % DPETSTEP MATLAB code for dPETSTEP.fig
    %      DPETSTEP, by itself, creates a new DPETSTEP or raises the existing
    %      singleton*.
    %
    %      H = DPETSTEP returns the handle to a new DPETSTEP or the handle to
    %      the existing singleton*.
    %
    %      DPETSTEP('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in DPETSTEP.M with the given input arguments.
    %
    %      DPETSTEP('Property','Value',...) creates a new DPETSTEP or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before dPETSTEP_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to dPETSTEP_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help dPETSTEP

    % Last Modified by GUIDE v2.5 10-Jul-2017 09:41:59

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @dPETSTEP_OpeningFcn, ...
                       'gui_OutputFcn',  @dPETSTEP_OutputFcn, ...
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

function dPETSTEP_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to dPETSTEP (see VARARGIN)

    % Choose default command line output for dPETSTEP
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

end

function varargout = dPETSTEP_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

function buttonSimulation_Callback(hObject, eventdata, handles)
    % hObject    handle to buttonSimulation (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    %% Run simulation GUI
    dPETSTEPgui_sim;

    %% Update handles structure
    guidata(hObject, handles);
end


function buttonModelFitting_Callback(hObject, eventdata, handles)
    % hObject    handle to buttonModelFitting (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    %% Run simulation GUI
    dPETSTEPgui_fit;

    %% Update handles structure
    guidata(hObject, handles);
end
