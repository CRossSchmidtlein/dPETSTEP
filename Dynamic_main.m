function [data,simSet,FBP4D,OS4D,OSpsf4D,counts,countsNoise] = Dynamic_main(data,frame,Cif,Cmean)
%******************************************************************************************************
% Run a complete dynamic PET simulation.
%
% USAGE  : [data,simSet,FBP4D,OS4D,OSpsf4D,counts,countsNoise] = Dynamic_main(data,frame,Cif,Cmean)
%
% INPUT  : data     Structure input data.
%          frame    Vector with start and end frame times in sec,
%                   [frameStart1; frameStart2=frameEnd1; frameStart3=frameEnd2;...].
%          Cif      Vector with input function to model (AIF or reference tissue TAC).
%          Cmean    Average activity for each frame of the simulation in (Bq/cc). 
%                   (Scales sinograms which in turn determines noise level).
%
% OUTPUT : data         Structure with all simulation input simulation data.
%          simSet       Structure with simulation settings.
%          FBP4D        Reconstructed dynamic FBP image in (Bq/cc). 
%          OS4D         Reconstructed dynamic OSEM image in (Bq/cc). 
%          OSpsf4D      Reconstructed dynamic OSEM w/ PSF image in (Bq/cc). 
%          counts       Pristine sinogram counts. 
%          countsNoise  Noisy sinogram counts. 
%******************************************************************************************************

%% Save log file.
logFile = 'Dynamic_main.log';
if exist(logFile, 'file')==2
  delete(logFile);
end
diary(logFile);
fprintf('%d-%02d-%02d, %02d:%02d:%02.0f\n',clock)

%% Timing.
mainClock = tic;

%% Sim settings.
simSet           = Dynamic_setSimParameters(frame,Cif,Cmean); %Initialize settings according to user specification

PIMscanNum       = simSet.PIMscanNum;      % PIM scan's ID.
PTscanNum        = simSet.PTscanNum;       % PT  scan's ID.
noFrames         = numel(frame)-1;         % Number of time points.
model			 = simSet.kineticModel;    % Desired kinetic model.
CifScaleFactor   = simSet.CifScaleFactor;  % Input function Cif scale factor. Multiply Cif with factor.
halflife         = simSet.halflife;        % Halflife of nuclide in sec. Can also be 'none' for no decay.
dt               = simSet.timeStep;        % Convolution time step in sec.
fovSize          = simSet.fovSize;         % Size of FOV. Voxel size = fovSize/simSize.
addVariability   = simSet.addVariability;  % Flag to add biologic variability (gaussian noise) or not.
variabilityScale = simSet.variabilityScale;% Scale factor of biologic variability.
FBP_OUT          = simSet.FBP_OUT;         % Flag to do FBP or not.
OS_OUT           = simSet.OS_OUT;          % Flag to do OSEM or not.
OSpsf_OUT        = simSet.OSpsf_OUT;       % Flag to do OSEM with PSF or not.

%% Mid frame time points.
midFrame         = frame(1:end-1) + diff(frame)/2;

%% Extract parametric image (PIM).
% Extract all parameter images (K1, k2,...) from 'data' and store in 4D PIM image
% where 4th dimension = parameter. Unit 1/sec for rate constants.
pim = data(2).data;
 
%% Add biologic variability to PIM (gaussian noise).
if addVariability
    pim = addVariability2PIM('pim',pim,'scale',variabilityScale);
end

%% Create pristine 4D PET image from parametric image (PIM).
% Arguments are the PIM, desired kinetic model, vector of frame start and
% stop times and arterial input function. Also, specify if wanting to use Matlabs parallell computing. 
littleClock = tic;
image4D     = createDynamicPETfromParametricImage_matrix('paramImage',pim,'model',model,...
    'frame',frame,'dt',dt,'Cif',Cif,'CifScaling',CifScaleFactor,...
    'doParallell',0,'doDecay',halflife); %Bq/cc

fprintf('\nTime for dynamic image generation: %.2f sec\n',toc(littleClock))

%% Pad 4D data with zeros or crop to get square and wanted recon voxel size and FOV.
voxSize = [ data(PIMscanNum).dataInfo.grid2Units data(PIMscanNum).dataInfo.grid1Units]; %unit (mm)
currFOV = [ voxSize(1)*size(image4D,1) voxSize(2)*size(image4D,2)]; %unit (mm)
image4D = padarray( image4D, [round((fovSize-currFOV(1))/2),round((fovSize-currFOV(2))/2),0,0] );

%% Add pristine 4D image to data structure
tmpStruct = data(PIMscanNum).dataInfo;
tmpStruct = rmfield(tmpStruct,'sizeOfDimension4');
data(PTscanNum).type                         = 'PT pristine';
data(PTscanNum).data                         = image4D; %Bq/cc
data(PTscanNum).Cif                          = Cif*CifScaleFactor; %Bq/cc
data(PTscanNum).midFrame                     = midFrame;
data(PTscanNum).dataInfo                     = tmpStruct;
data(PTscanNum).dataInfo.sizeOfDimension1    = size(image4D,1);
data(PTscanNum).dataInfo.sizeOfDimension2    = size(image4D,2);
data(PTscanNum).dataInfo.imageUnit           = 'BQML';
data(PTscanNum).dataInfo.numberOfDimensions  = 4;
data(PTscanNum).dataInfo.timeUnit            = 's';
data(PTscanNum).dataInfo.sizeOfDimensionTime = noFrames;
clear tmpStruct

%% Calculate voxel sizes, post filter, blurring kernels, mu-map etc.
% Doesn't run any simulations, just sets the joint values of simX (most
% values are the same for all frames).
[vox,PSFsim,PSFout,POST,scatterK,FWAC,initPT,sensScale] = Dynamic_PETSTEP_overhead(data,simSet);

%% Initialize.
if FBP_OUT
    FBP4D   = zeros( simSet.simSize, simSet.simSize, data(PTscanNum).dataInfo.sizeOfDimension3, noFrames, simSet.nREP ); 
else
    FBP4D   = [];
end
if OS_OUT    
    OS4D    = zeros( simSet.simSize, simSet.simSize, data(PTscanNum).dataInfo.sizeOfDimension3, noFrames, simSet.iterNUM, simSet.nREP ); 
else 
    OS4D    = [];
end
if OSpsf_OUT
    OSpsf4D = zeros( simSet.simSize, simSet.simSize, data(PTscanNum).dataInfo.sizeOfDimension3, noFrames, simSet.iterNUM, simSet.nREP ); 
else 
    OSpsf4D = [];
end
counts      = struct;
countsNoise = struct;
for i = 1:noFrames
    counts(i).total    = 0;
    counts(i).true     = 0;
    counts(i).scatter  = 0;
    counts(i).randoms  = 0;
    counts(i).NEC      = 0;
    counts(i).ID       = 0;
    countsNoise(i).total    = 0;
    countsNoise(i).true     = 0;
    countsNoise(i).scatter  = 0;
    countsNoise(i).randoms  = 0;
    countsNoise(i).NEC      = 0;
    countsNoise(i).ID       = 0;
end

%% Loop over all frames of the pristine 4D image and simulate each frame with PETSTEP.
littleClock = tic;
parfor i = 1:noFrames
    fprintf('\nFrame no %d/%d...\n',i,noFrames);
    
    output  = Dynamic_PETSTEP( data,simSet,i,vox,PSFsim,PSFout,POST,scatterK,FWAC,initPT,sensScale );
    
    % Assign output
    ind = 1;
    if FBP_OUT;   FBP4D(:,:,:,i,:)     = output{ind}; ind=ind+1; end
    if OS_OUT;    OS4D(:,:,:,i,:,:)    = output{ind}; ind=ind+1; end
    if OSpsf_OUT; OSpsf4D(:,:,:,i,:,:) = output{ind}; ind=ind+1; end
    counts(i)                          = output{ind}; ind=ind+1;
    countsNoise(i)                     = output{ind};
end
fprintf('\nTime for frame loop: %.2f min\n',toc(littleClock)/60)
clear output

%% Decay correct images
if ~strcmp(halflife,'none')
    fprintf('Decay correct reconstructed images...\n')
    lambda      = log(2)/halflife; %1/sec
    decayCorr   = (lambda*diff(frame)) ./ ( exp(-lambda*frame(1:end-1)) - exp(-lambda*frame(2:end)) );
    for i = 1:noFrames
        if FBP_OUT;   FBP4D(:,:,:,i,:)     = FBP4D(:,:,:,i,:)*decayCorr(i); end
        if OS_OUT;    OS4D(:,:,:,i,:,:)    = OS4D(:,:,:,i,:,:)*decayCorr(i); end
        if OSpsf_OUT; OSpsf4D(:,:,:,i,:,:) = OSpsf4D(:,:,:,i,:,:)*decayCorr(i); end
    end
end

%% End timing.
fprintf('\nTotal time: %.2f min\n',toc(mainClock)/60)

%% Clear variables.
% clear PTscanNum PIMscanNum model CifScaleFactor halflife fovSize dt addVariability variabilityScale FBP_OUT OS_OUT OSpsf_OUT 
% clear i logFile noFrames mainClock littleClock currFOV voxSize meanAct ind decayCorr doDecay lambda 
% clear vox PSFsim PSFout POST scatterK FWAC initPT sensScale
% clear pim image4D data

fprintf('%d-%02d-%02d, %02d:%02d:%02.0f\n',clock)
diary off

end