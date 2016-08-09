function output = Dynamic_PETSTEP(data,simSet,frameNo,vox,PSFsim,PSFout,POST,scatterK,FWAC,initPT,sensScale)
%%
%*******************************************************************************************************
%| Simulates an FBP, OSEM and OSEMpsf image from a pristine image.                                     |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| USAGE  :  [FBP,OS,OSpsf,counts,countsNoise] = Dynamic_PETSTEP(data,simSet,frameNo,vox,...           |
%|                             PSFsim,PSFout,POST,scatterK,FWAC,initPT)                                |
%|                                                                                                     |
%| INPUT  :  data           Structure with 4D PT and 3D CT data.                                       |
%|           simSet         Structure with all simulation settings.                                    |
%|           frameNo        Scalar with what frame from data structure to use.                         |
%|           vox            Structure with voxel information for input, simulation, and ouput.         |
%|           PSFsim         Matrix with PSF kernel matched to simulation size.                         |
%|           PSFout         Matrix with PSF kernel matched to output size.                             |
%|           POST           Matrix with postfilter kernel matched to output size.                      |
%|           scatterK       Matrix with scatter kernel.                                                |
%|           FWAC           Matrix with forward projection of mu-map.                                  |
%|           initPT         Matrix with initial PT guess for iterative recon.                          |
%|                                                                                                     |
%| OUTPUT :  FBP            Matrix with simulated FBP image.                                           |
%|           OS             Matrix with simulated OSEM image.                                          |
%|           OSpsf          Matrix with simulated PSF corrected OSEM image.                            |
%|           counts         Structure with pristine sinogram counts.                                   |
%|           countsNoise    Structure with noised sinogram counts.                                     |
%|                                                                                                     |
%*******************************************************************************************************
% Copyright 2016, C. Ross Schmidtlein, on behalf of the dPETSTEP development team.
% 
% This file is part of the Dynamic PET Simulator for Tracers via Emission Projection (dPETSTEP) software.
% 
% dPETSTEP development has been led by:  Ida H�ggstr�m, Bradley J. Beattie and C. Ross Schmidtlein.
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

%% input parameters
PTscanNum  = simSet.PTscanNum;       % PET scan's ID : should be automated

% Count Data
countScale = simSet.dwellTime(frameNo)*simSet.activityConc(frameNo)*simSet.countSens; % Sets mean number of counts per active voxel
SF         = simSet.SF;           % scatter fraction S/(T+S)
RF         = simSet.RF;           % randoms fraction R/(T+S+R)

% image reconstruction definitions
simSize    = simSet.simSize;         % Matrix size of reconstructed image
zFilter    = simSet.zFilter;         % post recon Z-axis filter 3-point smoothing
zFilter    = zFilter/sum(zFilter);   % automatic nomalization       
iterNUM    = simSet.iterNUM;         % number of iterations
subNUM     = simSet.subNUM;          % number of subsets

FBP_OUT    = simSet.FBP_OUT;
OS_OUT     = simSet.OS_OUT;
OSpsf_OUT  = simSet.OSpsf_OUT;

% number of replicate data sets
nREP       = simSet.nREP;

%% Get 3D PET and CT (might be mu-map) data
% Get 3D PET data (single frame)
PT = double(data(PTscanNum).data(:,:,:,frameNo));

%% Pad with zeros to square
sizePT     = size(PT);
maxDimPT   = max(sizePT(1:2));
padMatrix  = round( [maxDimPT - sizePT(1) maxDimPT - sizePT(2) 0] )/2;
PT         = padarray(PT,padMatrix);

%% Build simulated projection data
% build pristine data with count scaling
[FWPTtrue,FWPTscatter,FWPTrandoms,CTAC,wcc,counts] = Dynamic_buildSimFullData(...
    PT,FWAC,PSFsim,scatterK,vox,countScale,SF,RF,sensScale);

%% Allocate
if FBP_OUT
    FBP      = zeros([vox.petOut.nxn nREP]);
end
if OS_OUT
    OS       = zeros([vox.petOut.nxn iterNUM nREP]);
end
if OSpsf_OUT
    OSpsf    = zeros([vox.petOut.nxn iterNUM nREP]);
end

%% Perform Recons
for i = 1:nREP
    % add noise to projection data
    [nFWPTtotal,~,~,~,countsNoise] = ...
        noiseProjData(FWPTtrue,FWPTscatter,FWPTrandoms);
    countsNoise.ID      = countsNoise.NEC / prod( vox.petOut.nxn );

    % Recon Data
    if (FBP_OUT)
        [FBP(:,:,:,i)] = ...
            recon_fbpSimData_nonUniformSliceSens( ...
            nFWPTtotal,FWPTscatter+FWPTrandoms,CTAC, ...
            wcc,vox);
        [FBP(:,:,:,i)] = xyPostFilter(FBP(:,:,:,i),POST);
        [FBP(:,:,:,i)] = zAxialFilter(FBP(:,:,:,i),zFilter);
    end
    if (OS_OUT)
        [OS(:,:,:,:,i)] = ...
            recon_osemSimData_nonUniformSliceSens( ...
            nFWPTtotal,FWPTscatter+FWPTrandoms,CTAC,initPT, ...
            wcc,vox,iterNUM,subNUM);
        [OS(:,:,:,:,i)] = xyPostFilter(OS(:,:,:,:,i),POST);
        [OS(:,:,:,:,i)] = zAxialFilter(OS(:,:,:,:,i),zFilter);
    end
    if (OSpsf_OUT)
        [OSpsf(:,:,:,:,i)] = ...
            recon_osemPSFSimData_nonUniformSliceSens( ...
            nFWPTtotal,FWPTscatter+FWPTrandoms,CTAC,initPT,PSFout, ...
            wcc,vox,iterNUM,subNUM);
        [OSpsf(:,:,:,:,i)] = xyPostFilter(OSpsf(:,:,:,:,i),POST);
        [OSpsf(:,:,:,:,i)] = zAxialFilter(OSpsf(:,:,:,:,i),zFilter);
    end

end

% Output: Recons, counts and wcc
ind = 1;
if FBP_OUT
    output{ind} = FBP;
    ind = ind + 1;
end
if OS_OUT
    output{ind} = OS;
    ind = ind + 1;
end
if OSpsf_OUT
    output{ind} = OSpsf;
    ind = ind + 1;
end
output{ind}   = counts;      ind = ind + 1;
output{ind}   = countsNoise; ind = ind + 1;
output{ind}   = nFWPTtotal;  ind = ind + 1;
output{ind}   = FWPTtrue;    ind = ind + 1;
output{ind}   = FWPTscatter; ind = ind + 1;
output{ind}   = FWPTrandoms; ind = ind + 1;
output{ind}   = wcc;

return
