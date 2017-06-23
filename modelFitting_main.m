function paramImage = modelFitting_main(varargin)
%%
%****************************************************************************************************
%| Fits the dynamic PET data to a kinetic model, producing parametric images (voxel-wise) or a set  |
%| of parameters for a ROI.                                                                         |
%|                                                                                                  |
%| IH, 19/04/2016                                                                                   |
%|                                                                                                  |
%| N = no of kinetic parameters.                                                                    |
%| f = no of frames.                                                                                |
%|                                                                                                  |
%| USAGE  :   paramImage = modelFitting_main('image',image,...                                      |
%|                                    'model','2tissue',...                                         |
%|                                    'midFrame',[f_1 f_2...f_f-1],...                              |
%|                                    'w',[w_1 w_2...w_f-1] or e.g. 'w1', 'ones',...                |
%|                                    'p0',[p0_1 p0_2 ... p0_N],...                                 |
%|                                    'ROIMask',[nx,ny,nz],...                                      |
%|                                    'Cp',[cp_1 cp_2...cp_f-1],...                                 |
%|                                    'CpMask',[nx,ny,nz],...                                       |
%|                                    'Cref',[c_1 c_2...c_f-1],...                                  |
%|                                    'CrefMask',[nx,ny,nz],...                                     |
%|                                    'halflife', scalar in sec,...                                 |
%|                                    'noExp',scalar,...                                            |
%|                                    'interpMethod',string,...                                     |
%|                                    'solver',string,...                                           |
%|                                    'lowerBound', [lb_1 lb_2...lb_N],...                          |
%|                                    'upperBound', [ub_1 ub_2...ub_N],...                          |
%|                                    'noCPU', noCPU)                                               |
%|                                                                                                  |
%| INPUT  :   image         4D matrix with dynamic image, unit (arbitrary), e.g. (Bq/ml).           |
%|                          image = [nx,ny,nz,f-1].                                                 |
%|            model         String of desired kinetic model, '1Tissue', '2Tissue', 'FRTM', 'SRTM',  |
%|                          or 'sumExp'.                                                            |
%|            midFrame      Vector with mid frame times, unit (s).                                  |
%|                          midFrame = [f-1,1].                                                     |
%|            Cp            Vector with blood input function to model, unit (arbitrary), e.g.       |
%|                          (Bq/cc).                                                                |
%|                          Cp = [(f-1),1].                                                         |
%|            Cref          Vector with reference tissue TAC, unit (arbitrary), e.g. (Bq/cc). Used  |
%|                          for kinetic models FRTM and SRTM.                                       |
%|                          Cref = [(f-1),1].                                                       |
%|                          Cp = [(f-1),1].                                                         |
%|            OPTIONAL :                                                                            |
%|            ROIMask       3D matrix with image ROI mask, for an image-derived response function.  |
%|                          ROIMask = [nx,ny,nz].                                                   |
%|            CpMask        3D matrix with Cp ROI mask, for an image-derived input function.        |
%|                          CpMask = [nx,ny,nz].                                                    |
%|            CrefMask      3D matrix with Cref ROI mask, for an image-derived Cref.                |
%|                          CrefMask = [nx,ny,nz].                                                  |
%|            w             String or vector with frame weights, unit (unitless).                   |
%|                          w = [(f-1),1], or string 'w1', 'w2',...,'w6' or 'ones' (default).       |
%|            p0            Vector with initial parameter guesses, unit (1/s) for rate constants.   |
%|                          p0 = [1,N]. Default 0.01 for all parameters.                            |
%|            halflife      Nuclide halflife in sec, for calculation of some weight types.          |
%|            noExp         Scalar between 1-inf. Number of exponentials for model 'sumExp'.        |
%|            interpMethod  String with interpolation method. Default 'linear'.                     |
%|            solver        String with desired solver algorith. Default []-->'trust-region-refl'.  |
%|            lowerBound    Vector with lower bound for parameters, unit (1/s) for rate constants.  |
%|                          lowerBound = [1,N]. Default 0 for all parameters.                       |
%|            upperBound    Vector with upper bound for parameters, unit (1/s) for rate constants.  |
%|                          upperBound = [1,N]. Default 100*p0.                                     |
%|            noCPU         Number of CPUs to use (for voxelwise fitting). Default 1.               |
%|                                                                                                  |
%| OUTPUT :   paramImage    1) No ROIMask: 4D matrix with parametric image (3D image with a 4th     |
%|                          parameter dimension), unit (1/s) for rate constants.                    |
%|                          paramImage = [nx,ny,nz,N].                                              |
%|                          2) ROIMask specified: vector of parameters for ROI.                     |
%|                          paramImage = [1,N].                                                     |
%|                                                                                                  |     
%****************************************************************************************************
% Copyright 2016, C. Ross Schmidtlein, on behalf of the dPETSTEP development team.
% 
% This file is part of the Dynamic PET Simulator for Tracers via Emission Projection (dPETSTEP) software.
% 
% dPETSTEP development has been led by:  Ida Häggström, Bradley J. Beattie and C. Ross Schmidtlein.
% 
% dPETSTEP has been financially supported by the Cancer Research Foundation in Northern Sweden, 
% the US National Institutes of Health and the National Cancer Institute under multiple grants.
% 
% dPETSTEP is distributed under the terms of the Lesser GNU Public License. 
% 
%     This version of dPETSTEP is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
% dPETSTEP is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with dPETSTEP.  If not, see <http://www.gnu.org/licenses/>.
%******************************************************************************************************

%% Read all function arguments.
for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 'image'
            image     = varargin{i+1};  % 4D matrix with dynamic image.
                                        % image = [nx,ny,nz,f-1], unit (arbitrary) e.g. (Bq/ml).
        case 'midFrame'
            midFrame  = varargin{i+1};  % Mid frame times in (sec). 
        case 'Cp'
            Cp        = varargin{i+1};  % Input function to the kinetic model, in unit (Bq/ml). 
                                        % Cp = [(f-1),1], f = no of frames.
        case 'Cref'
            Cref      = varargin{i+1};  % Reference tissue TAC, in unit (Bq/ml). 
                                        % Cref = [(f-1),1], f = no of frames.
        case 'ROIMask'
            ROIMask   = varargin{i+1};  % 3D matrix with image ROI mask, for an image-derived response function. 
                                        % ROIMask = [nx,ny,nz].
        case 'CpMask'
            CpMask    = varargin{i+1};  % 3D matrix with Cp ROI mask, for an image-derived Cp.  
                                        % CpMask = [nx,ny,nz].
        case 'CrefMask'
            CrefMask  = varargin{i+1};  % 3D matrix with Cref ROI mask, for an image-derived Cref.  
                                        % CrefMask = [nx,ny,nz].
        case 'model'
            modelName = varargin{i+1};  % Kinetic model, e.g. '2-tissue'.
        case 'p0'                          
            p0        = varargin{i+1};  % Initial guesses for parameters
                                        % p0 = [1,noP], noP = no of parameters.
        case 'w'
            wIn       = varargin{i+1};  % Vector with frame weights.
                                        % w = [(f-1),1], f = no of frames.
        case 'halflife'
            halflife  = varargin{i+1};  % Scalar, nuclide halflife in sec.
        case 'noExp'
            noExp     = varargin{i+1};  % Scalar (1-inf), number of exponentials for model "sumExp".
        case 'interpMethod'
            interpMethod = varargin{i+1};  % String with interpolation method.
        case 'solver'
            solver    = varargin{i+1};  % String with desired solver algorithm. Defaults to []-->"trust-region-reflective".
        case 'lowerBound'
            lowerBound   = varargin{i+1};  % Vector with lower bounds for parameters.
        case 'upperBound'
            upperBound   = varargin{i+1};  % Vector with upper bounds for parameters.
        case 'noCPU'
            noCPU   = varargin{i+1};  % Scalar with number of CPUs to use (Default 1).
        otherwise
            fprintf('Unknown argument ''%s''.\nExiting...\n',varargin{i});
            paramImage = []; 
            return;
    end
end

%% No initial guess.
if ~exist('p0','var')
    fprintf('Start guess (p0) required to know number of parameters!\n')
    paramImage = [];
    return;
end

%% Number of kinetic parameters.
noP           = numel(p0);

%% Default upper and lower bounds for fit.
if ~exist('lowerBound','var')
    fprintf('DEFAULT: Setting lower parameter bounds to all zeros...\n')
    lowerBound = zeros(size(p0));
end
if ~exist('upperBound','var')
    fprintf('DEFAULT: Setting upper parameter bounds 100*p0...\n')
    upperBound = 100*p0;
end
if ~exist('solver','var')
    solver = [];
end

%% No of CPUs.
if ~exist('noCPU','var')
    fprintf('DEFAULT: Using single CPU.\n')
    parforArg = 0;
elseif nCPU==1
    fprintf('Using single CPU.\n')
    parforArg = 0;
else
    fprintf('Using %d CPUs.\n',noCPU)
    parforArg = noCPU;
end

%% Default interpolation method.
if ~exist('interpMethod','var')
    fprintf('DEFAULT: Setting interpolation method to ''linear''...\n')
    interpMethod = 'linear'; % default
end

%% Threshold level. Skip fit for activity below this. User can change level here.
threshold      = 0.01; %unitless. 

%% Dynamic image size.
sizeIm         = size(image);
noVoxels       = sizeIm(1)*sizeIm(2)*sizeIm(3);

%% Number of exponentials if model is "sumExp"
if strcmp('model','sumExp') && ~exist('noExp','var')
    fprtinf('You need to specify "noExp" for model "sumExp"!\n')
    paramImage = [];
    return;
end

%% ROI mask for Cp or Cref.
if exist('CpMask','var')
    fprintf('Using ROI mask for image-derived Cp...\n')
    Cp   = pickOutTACFromROI(image,CpMask,'average');
end
if exist('CrefMask','var')
    fprintf('Using ROI mask for image-derived Cref...\n')
    Cref = pickOutTACFromROI(image,CrefMask,'average');
end

%% Calculate or use existing frame weights.
if ~exist('wIn','var')
    fprintf('DEFAULT: Setting frame weights to uniformly 1 for all frames...\n')
    w            = calculateWeights('t',midFrame,'type','ones');
elseif ischar(wIn)
    fprintf('Calculating frame weights according to "%s"...\n',wIn)
    if exist('halflife','var')
        w        = calculateWeights('t',midFrame,'type',wIn,'halflife',halflife);
    else
        w        = calculateWeights('t',midFrame,'type',wIn);
    end
elseif isvector(wIn)
    fprintf('Using supplied frame weights...\n')
    w            = wIn;
end

%% Settings for chosen model.
switch modelName
    case {'1Tissue','1tissue','1-tissue','1-Tissue'}
        inputFunc = 'Cp'; %model input function
        fitFunc   = @(t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm)   fit_1Tissue_lsqnonlin( t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm); %fit funtion
    case {'2Tissue','2tissue','2-tissue','2-Tissue'}
        inputFunc = 'Cp';
        fitFunc   = @(t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm)   fit_2Tissue_lsqnonlin( t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm); %fit funtion
    case 'FRTM'
        inputFunc = 'Cref';
        fitFunc   = @(t,C,Cref,w,p0,dt,doPlot,LB,UB,algorithm) fit_FRTM_lsqnonlin(    t,C,Cref,w,p0,dt,doPlot,LB,UB,algorithm); %fit funtion
    case 'SRTM'
        inputFunc = 'Cref';
        fitFunc   = @(t,C,Cref,w,p0,dt,doPlot,LB,UB,algorithm) fit_SRTM_lsqnonlin(    t,C,Cref,w,p0,dt,doPlot,LB,UB,algorithm); %fit funtion
    case {'sumExp','sumexp'}
        inputFunc = 'Cp';
        fitFunc   = @(t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm)   fit_sumExp_lsqnonlin(  t,C,Cp,w,p0,dt,doPlot,LB,UB,algorithm); %fit funtion
end

%% Interpolate to equidistant time step dt (for convolution).
dt                 = min(diff(midFrame));
t2                 = (midFrame(1):dt:midFrame(end))';
w2                 = interpolateWeights(midFrame,w,t2);
if exist('Cp','var')
    Cp2            = interp1(midFrame,Cp,t2,interpMethod,'extrap');
    Cp2(Cp2<0)     = 0;
    inputFunc      = Cp2;
end
if exist('Cref','var')
    Cref2          = interp1(midFrame,Cref,t2,interpMethod,'extrap');
    Cref2(Cref2<0) = 0;
    inputFunc      = Cref2;
end
    
%% Fit to model, ROI-wise or voxel-wise.
if exist('ROIMask','var') % ROI-WISE ######################################
    fprintf('Doing single ROI fitting...\n')  
    C           = pickOutTACFromROI(image,ROIMask,'average');
    C2          = interp1(midFrame,C,t2,interpMethod,'extrap');
    C2(C2<0)    = 0;
    paramImage  = fitFunc(t2,C2,inputFunc,w2,p0,dt,1,lowerBound,upperBound,solver);
else %VOXEL-WISE ##########################################################
    fprintf('Doing voxel-wise fitting of %d voxels...\n',noVoxels)   
    %% Total activity of Cp. Used for thresholding.
    act         = sum(Cp);   
    %% Reshape 4D image to 2D matrix [noVoxels,noTimePoints].
    image       = reshape( image, [noVoxels,numel(midFrame)]);
    %% Allocate parameter matrix.
    p           = zeros(noVoxels,noP);
    counter     = zeros(noVoxels,1);
    %% Loop over all voxels.
    tic
    parfor (k = 1:noVoxels, parforArg)
        if sum(image(k,:))/act > threshold 
            counter(k) = 1;
            C2         = interp1(midFrame,image(k,:)',t2,interpMethod,'extrap');
            C2(C2<0)   = 0;
            p(k,:)     = fitFunc(t2,C2,inputFunc,w2,p0,dt,0,lowerBound,upperBound,solver);
        end
    end
    fprintf('%d voxels of %d above threshold (%.0f%%).\n',sum(counter), noVoxels, 100*sum(counter)/noVoxels)
    fprintf('Loop time : %.2f minutes (%.2f seconds).\n',toc/60, toc)
    
    %% Reshape back to image dimensions.
    paramImage    = reshape(p,[sizeIm(1) sizeIm(2) sizeIm(3) noP]);
end
