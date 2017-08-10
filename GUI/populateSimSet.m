function simSet = populateSimSet(handles)

%% Halflife
if isempty(handles.halflife.String) || strcmp(handles.halflife.Enable,'off')
    halflife = 'none';
else
    halflife = handles.halflife.Value;
end

%% Interp method
names = get(findobj(allchild(handles.Interp),'style','radiobutton'),'Tag');
for i=1:numel(names)
    pushed = get(handles.(names{i}),'value');
    if pushed
        interpMethod = handles.(names{i}).String;
    end
end

%% Kinetic model
names = get(findobj(allchild(handles.KineticModel),'style','radiobutton'),'Tag');
for i=1:numel(names)
    pushed = get(handles.(names{i}),'value');
    if pushed
        model = handles.(names{i}).String;
    end
end

%% Number of rad bin directly or from scanner diameter
if handles.fromRadBin.Value==1 %directly
    numRadBin = handles.numRadBin.Value;
    ringData  = [];
else % from scanner diameter
    numRadBin = [];
    ringData  = handles.ringDiameter.Value;
end

%% Read frame from WS or file
if handles.frameFromWS.Value==1 %selected read from WS
    frame   = evalin('base',handles.frameName.String);
else
    [~,~,ext] = fileparts(handles.frameName.String);
    if strcmp(ext,'.mat') %MAT-file
        frame = cell2mat(struct2cell( load(handles.frameName.String) ));
    else %Text file
        frame = load(handles.frameName.String);
    end
end
if numel(frame)<1
    fprintf('File "%s" is empty.\n',handles.frameName.String);
    simSet = struct;
    return
elseif numel(frame)>1
    dwellTime = diff(frame);
else
    dwellTime = 1;
end

%% Read input function from WS or file
if handles.ifFromWS.Value==1 %selected read from WS
    Cif   = evalin('base',handles.inputFuncName.String);
else
    [~,~,ext] = fileparts(handles.inputFuncName.String);
    if strcmp(ext,'.mat') %MAT-file
        Cif = cell2mat(struct2cell( load(handles.inputFuncName.String) ));
    else %Text file
        Cif = load(handles.inputFuncName.String);
    end
end

%% Z-postfilter
if strcmp(handles.postfilterZ.Enable,'off')
    zFilter = [];
else
    if handles.zfilterFromWS.Value==1 %selected read from WS
        zFilter = evalin('base',handles.postfilterZ.String);
    else
        [~,~,ext] = fileparts(handles.postfilterZ.String);
        if strcmp(ext,'.mat') %MAT-file
            zFilter = cell2mat(struct2cell( load(handles.postfilterZ.String) ));
        else %Text file
            zFilter = load(handles.postfilterZ.String);
        end
    end
end

%% Simulation parameters
simSet = struct;

% count data
simSet.countSens      = handles.countSens.Value; % Scanner sensitivity (counts/kBq/s)
simSet.SF             = handles.SF.Value;        % scatter fraction S/(T+S)
simSet.RF             = handles.RF.Value;        % randoms fraction R/(T+S+R)

% Dynamic settings
simSet.frame          = frame;     % Vector with frame times in sec.
simSet.dwellTime      = dwellTime; % Frame lengths.
simSet.Cif            = Cif;       % Vector with input function in Bq/cc. Either arterial input function.
                                   % or reference tissue TAC, depending on what model you use.
simSet.CifScaleFactor = str2num(handles.CifScaleFactor.String); % Scale factor to multiply the supplied input function with.
simSet.halflife       = halflife;                               % Halflife of nuclide in sec, or 'none' for no decay.
simSet.timeStep       = str2num(handles.timeStep.String);       % Convolution time step in sec.
simSet.interpMethod   = interpMethod;                           % Interpolation method.
simSet.kineticModel   = model;                                  % Desired kinetic model, '1-Tissue', '2-Tissue', 'FRTM', 'SRTM' or 'sumExp'.

% scanner charaterisitics
simSet.RingData    = ringData;                          % the diameter of the scanner ring
simSet.radBin      = numRadBin;                         % Sets inital radial projetion data size
simSet.tanBin      = handles.numAngBin.Value;           % Sets inital angular projetion data size
simSet.maxRingDiff = handles.maxRingDiff.Value;         % Maximum allowed ring difference
simSet.psf         = str2num(handles.psf.String);       % Correction PSF.
simSet.blurT       = handles.blurT.Value;               % System blurring PSF

% image reconstruction definitions
simSet.fovSize    = handles.fovSize.Value;       % size of dynamic image FOV in mm. Voxel size = fovSize/simSize.
simSet.simSize    = handles.simSize.Value;       % matrix size of reconstructed image
simSet.zFilter    = zFilter;                     % post recon Z-axis filter 3-point smoothing
                                                 % Heavy[ 1 2 1]/4, Standard[1 4 1]/6, Light[1 6 1]/8, None[0 1 0]
simSet.postFilter = [];
simSet.iterNUM    = [];
simSet.subNUM     = [];
if strcmp(handles.postfilterXY.Enable,'on')
    simSet.postFilter = str2num(handles.postfilterXY.String); % FWHM in (mm) of post reconstruction filter. FWHM = 2*sqrt(2*log(2))*sigma.
end
if strcmp(handles.iterNUM.Enable,'on')
    simSet.iterNUM    = str2num(handles.iterNUM.String);      % number of iterations
end
if strcmp(handles.subNUM.Enable,'on')
    simSet.subNUM     = str2num(handles.subNUM.String);       % number of subsets
end

% Biologic variability
simSet.addVariability   = handles.addBioVar.Value;           % Flag to add biologic variability or not
simSet.variabilityScale = str2num(handles.varScale.String);  % Scale factor of variability (variance of gaussian noise = image/scale)

% Reconstructions
simSet.FBP_OUT    = handles.doFBP.Value;      % Do FBP recon
simSet.OS_OUT     = handles.doOSEM.Value;     % Do OSEM recon
simSet.OSpsf_OUT  = handles.doOSEMPSF.Value;  % Do OSEM w/ PSF recon

% number of replicate data sets
simSet.nREP       = handles.noREP.Value;         % Number of noise realizations

end