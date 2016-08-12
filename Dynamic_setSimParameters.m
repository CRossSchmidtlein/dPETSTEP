function simSet = Dynamic_setSimParameters(frame,Cif)
%******************************************************************************************************
% Sets dPETSTEP simulation parameters.
%
% USAGE  : simSet = Dynamic_setSimParameters(frame,Cif)
%
% INPUT  : frame    Vector with start and end frame times in sec,
%                   [frameStart1; frameStart2=frameEnd1; frameStart3=frameEnd2;...].
%          Cif      Vector with input function to model (AIF or reference tissue TAC).
%
% OUTPUT : simSet  Structure with all simulation settings
%******************************************************************************************************
% IH, 19/04/2016
%
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

%% Simulation parameters
simSet = struct;
 
% input parameters
simSet.CTscanNum  = 1;                  % CT  scan's ID : should be automated
simSet.PIMscanNum = 2;                  % Parametric image scan's ID : should be automated
simSet.PTscanNum  = 3;                  % Dynamic PET image scan's ID : should be automated

% count data
simSet.countSens      = (265/324)*6.44; % 3D sensitivity (counts/kBq/s) (GE DLS)
simSet.SF             = 0.289;          % scatter fraction S/(T+S)
simSet.RF             = 0.02;           % randoms fraction R/(T+S+R)

% Dynamic settings
simSet.frame          = frame;          % Vector with frame times in sec.
simSet.dwellTime      = diff(frame);    % Frame lengths.
simSet.Cif            = Cif;            % Vector with input function in Bq/cc. Either arterial input function.
                                        % or reference tissue TAC, depending on what model you use.
simSet.CifScaleFactor = 1;              % Scale factor to multiply the supplied input function with.
simSet.halflife       = 'none';         % Halflife of nuclide in sec, or 'none' for no decay.
simSet.timeStep       = 0.5;            % Convolution time step in sec.
simSet.interpMethod   = 'linear';        % Interpolation method.
% simSet.kineticModel   = '2-Tissue';		% Desired kinetic model, '1-Tissue', '2-Tissue', 'FRTM', 'SRTM' or 'sumExp'.
simSet.kineticModel   = 'sumExp';		% Desired kinetic model, '1-Tissue', '2-Tissue', 'FRTM', 'SRTM' or 'sumExp'.

% scanner charaterisitics
simSet.RingData   = 430%880;        % the diameter of the scanner ring (GE DLS)
simSet.tanBin     = 336;        % Sets inital projetion data size (GE DLS)
simSet.maxRingDiff= 11;         % Maximum allowed ring difference
simSet.psf        = 5.1;        % Assumes a PSF for the system, uses same for correction. FWHM.
simSet.blurT      = 5.1;        % DLS

% image reconstruction definitions
simSet.fovSize    = 128;        % size of dynamic image FOV in mm. Voxel size = fovSize/simSize.
simSet.simSize    = 128;        % matrix size of reconstructed image
simSet.zFilter    = [1];    % post recon Z-axis filter 3-point smoothing
                                % Heavy[ 1 2 1]/4, Standard[1 4 1]/6, Light[1 6 1]/8, None[0 1 0]
simSet.postFilter = 6;          % FWHM in (mm) of post reconstruction filter. FWHM = 2*sqrt(2*log(2))*sigma.
simSet.iterNUM    = 5;          % number of iterations
simSet.subNUM     = 16;         % number of subsets

% Biologic variability
simSet.addVariability   = false;  % Flag to add biologic variability or not
simSet.variabilityScale = 10;     % Scale factor of variability (variance of gaussian noise = image/scale)

% Reconstructions
simSet.FBP_OUT       = 0;
simSet.OS_OUT        = false;
simSet.OSpsf_OUT     = false;

% number of replicate data sets
simSet.nREP       = 1;

end