function [vox,PSFsim,PSFout,POST,scatterK,FWAC,initPT,sensScale,activityConc] = Dynamic_PETSTEP_overhead(data,simSet)
%%
%******************************************************************************************************
% Calculates voxel sizes, attenuation factors, initial PET etc.
%
% USAGE  : [vox,PSFsim,PSFout,POST,scatterK,FWAC,initPT,sensScale] = Dynamic_PETSTEP_overhead(data,simSet)
%
% INPUT  : data      Structure with CT/mumap and dynamic image
%          simSet    Structure with all simulation settings
%
% OUTPUT : vox  	 structure with voxel sizes
%		   PSFsim	 Matrix with system PSF kernel
%		   PSFout	 Matrix with PSF kernel for correction during reconstruction
%          POST  	 Matrix with Gaussian XY post filter kernel
% 	       scatterK  Matrix with scatter kernel
%	       FWAC		 Matrix with attenuation factors in sinogram space
%          initPT  	 Initial guess of PET image (disk of ones)
%		   sensScale Vector with sensitivity scale factor for all slices
%******************************************************************************************************
% IH, 19/04/2016
%
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
CTscanNum  = simSet.CTscanNum;       % CT  scan's ID : should be automated
noFrames   = size(data(CTscanNum).data,4);

% scanner charaterisitics
ringData   = simSet.RingData;     % the ring diameter of the scanner's data acquisition (810 mm GE D690)
tanBin     = simSet.tanBin;       % Sets inital projetion data size (280 matched to GE DSTE)
psf        = simSet.psf;          % Assumes a PSF for the system, uses same for correction

% image reconstruction definitions
simSize    = simSet.simSize;         % Matrix size of reconstructed image 
postFilter = simSet.postFilter;      % fwhm of post reconstruction filter


% Convert CT to mu-map if not already mu-map
switch data(CTscanNum).type
    case 'CT'
        % Get CT (or mu-map) data
        CT              = double(data(CTscanNum).data);
        CT(CT < -500) = -1000;
        CTmuTmp         = convertCT2mumap(CT);
    case 'mumap'
        CTmuTmp         = double(data(CTscanNum).data);
end

%% Get inital image and data attributes 
% Get reference PET voxel size in mm and data
vox.pet.xyz  = [ data(PTscanNum).dataInfo.grid2Units ...
                data(PTscanNum).dataInfo.grid1Units ...
                data(PTscanNum).dataInfo.grid3Units ];
vox.pet.vol  = prod(vox.pet.xyz);
vox.pet.nxn  = [ data(PTscanNum).dataInfo.sizeOfDimension1 ...
                data(PTscanNum).dataInfo.sizeOfDimension2 ...
                data(PTscanNum).dataInfo.sizeOfDimension3 ];
vox.pet.fov  = vox.pet.nxn(1:3).*vox.pet.xyz;
vox.pet.rtz  = [ 2*ceil(norm( vox.pet.nxn(1:2) - floor(( vox.pet.nxn(1:2)-1 )/2)-1)) + 3 ...
                tanBin vox.pet.nxn(3) ];
vox.pet;

% Get reference CT voxel size in mm and data
vox.muct.xyz = [ data(CTscanNum).dataInfo.grid2Units ...
               data(CTscanNum).dataInfo.grid1Units ...
               data(CTscanNum).dataInfo.grid3Units ];
vox.muct.vol = prod(vox.muct.xyz);
vox.muct.nxn = [ data(CTscanNum).dataInfo.sizeOfDimension1 ...
                data(CTscanNum).dataInfo.sizeOfDimension2 ...
                data(CTscanNum).dataInfo.sizeOfDimension3 ];
% vox.muct.nxn = size(CTmuTmp);
vox.muct.fov = vox.muct.nxn.*vox.muct.xyz;

% Get Simulation PET/CT voxel size in mm and data
radBin = floor(4*tanBin * asin( vox.pet.fov(1) / ringData )/pi - 1);
if (mod(radBin,2) == 0), radBin = radBin + 1 ; end
vox.petSim.xyz = vox.pet.xyz .* [ vox.pet.nxn(1:2)/radBin 1 ];
vox.petSim.vol = prod(vox.petSim.xyz);
vox.petSim.nxn = [radBin radBin vox.pet.nxn(3)];
vox.petSim.fov = vox.petSim.nxn.*vox.petSim.xyz;
vox.petSim.r   = 2*ceil(norm( vox.petSim.nxn(1:2) - floor(( vox.petSim.nxn(1:2)-1 )/2)-1)) + 3;
vox.petSim.tan = tanBin;
vox.petSim.rtz = [ 2*ceil(norm( vox.petSim.nxn(1:2) - floor(( vox.petSim.nxn(1:2)-1 )/2)-1)) + 3 ...
                   tanBin vox.pet.nxn(3) ];
vox.petSim

% Get Output Image PET voxel size in mm and data
vox.petOut.xyz = vox.pet.xyz .* [ vox.pet.nxn(1:2)/simSize 1 ];
vox.petOut.vol = prod(vox.petOut.xyz);
vox.petOut.nxn = [simSize simSize vox.pet.nxn(3)];
vox.petOut.fov = vox.petOut.nxn.*vox.petOut.xyz;
vox.petOut.rtz = [ 2*ceil(norm( vox.petOut.nxn(1:2) - floor(( vox.petOut.nxn(1:2)-1 )/2)-1)) + 3 ...
                   tanBin vox.pet.nxn(3) ];
vox.petOut;

%% Z-axis slice sentivity correction
if (mod(vox.petSim.rtz(3),2) == 1)
    sensScale = [1:floor(vox.petSim.rtz(3)/2)+1 floor(vox.petSim.rtz(3)/2):-1:1]; 
else
    sensScale = [1:floor(vox.petSim.rtz(3)/2) floor(vox.petSim.rtz(3)/2):-1:1]; 
end
sensScale(sensScale>(simSet.maxRingDiff+1)) = max( sensScale( sensScale<=(simSet.maxRingDiff+1) ) );
sensScale = sensScale/mean(sensScale);

% PSF kernel matched to simulation size
fwhm        = psf * vox.petSim.nxn(1) / vox.pet.fov(1);
fwhmMat     = max( ceil(3*fwhm) , 5 );
if (mod(fwhmMat,2) == 0)
    fwhmMat = fwhmMat + 1; 
end
PSFsim     = fspecial('gaussian', fwhmMat, fwhm / ( 2*sqrt(2*log(2)) ) );

% PSF kernel matched to output size
fwhm        = psf * vox.petOut.nxn(1) / vox.pet.fov(1);
fwhmMat     = max( ceil(3*fwhm) , 5 );
if (mod(fwhmMat,2) == 0), fwhmMat = fwhmMat + 1; end
PSFout     = fspecial('gaussian', fwhmMat, fwhm / ( 2*sqrt(2*log(2)) ) );

% Post smoothing kernel matched to output size
fwhmPost    = postFilter * vox.petOut.nxn(1) / vox.pet.fov(1);
fwhmMatPost = max( ceil(3*fwhmPost) , 5 );
if (mod(fwhmMatPost,2) == 0), fwhmMatPost = fwhmMatPost + 1; end
POST        = fspecial('gaussian', fwhmMatPost, fwhmPost / ( 2*sqrt(2*log(2)) ) );

% scatter kernel
scatterFWHM  = 200; %(mm)
fwhmS        = scatterFWHM * vox.petSim.nxn(1) / vox.pet.fov(1);
fwhmMat      = max( ceil(3*fwhmS) , 3 );
if (mod(fwhmMat,2) == 0)
    fwhmMat  = fwhmMat + 1; 
end
scatterK     = fspecial('gaussian',fwhmMat, fwhmS / ( 2*sqrt(2*log(2)) ) );

%% Rescale mumap image
CTmu   = zeros([vox.petSim.nxn noFrames]);
for j = 1:noFrames
    for i = 1:vox.petSim.nxn(3)
        if (isempty(CTmuTmp))
            CTmu(:,:,:,j) = [];
        else
            % pad CT to match PET
            tmp0  = zeros(vox.petSim.nxn(1:2));
            tmp1  = imresize(CTmuTmp(:,:,i,j),vox.muct.xyz(1)/vox.petSim.xyz(1),'cubic');
            
            xA = max(1,round(( vox.petSim.nxn(1) - size(tmp1,2) )/2) + 1);
            xB = xA + size(tmp1,2) - 1;
            yA = max(1,round(( vox.petSim.nxn(1) - size(tmp1,1) )/2) + 1);
            yB = yA + size(tmp1,1) - 1;
            
            tmp0(yA:yB,xA:xB) = tmp1;
            CTmu(:,:,i,j) = imresize( tmp0,vox.petSim.nxn(1:2), 'cubic' );
        end
    end
end
CTmu(CTmu < 0)     = 0; % Alignment verified
% Blur with PSF
for j = 1:noFrames
    for i = 1:vox.petSim.nxn(3)
        CTmu(:,:,i,j)    = imfilter(CTmu(:,:,i,j),   PSFsim,      'replicate','same','conv');
    end
end

%% Forward project attenuation image (mumap)
PHI         = 0:180/vox.petSim.rtz(2):180*(1-1/vox.petSim.rtz(2));
FWAC        = ones([vox.petSim.rtz(1) vox.petSim.rtz(2) vox.petSim.rtz(3) noFrames]);
for j = 1:noFrames
    for i = 1:vox.petSim.rtz(3)
        if (isempty(CTmu))
            FWAC(:,:,i,j) = ones([vox.petSim.rtz(1) vox.petSim.rtz(2)]);
        else
            FWAC(:,:,i,j) = exp( -( vox.pet.xyz(1) * vox.pet.nxn(1) / vox.petSim.nxn(1) )/10 * ...
                radon( CTmu(:,:,i,j), PHI ) );
        end % Alignment verified
    end
end
FWAC(FWAC > 1) = 1;

%Initial PET
initPT = zeros( vox.petOut.nxn );
nn     = (vox.petOut.nxn(1)-1)/2;
x2     = linspace(-nn,nn,vox.petOut.nxn(1)).^2;
disk   = x2(ones(vox.petOut.nxn(1),1),:) + x2(ones(1,vox.petOut.nxn(1)),:)' <= (vox.petOut.nxn(1)/2)^2;
for i  = 1:vox.petOut.nxn(3)
    initPT(:,:,i) = disk;
end

% Attenuationed total activity per frame
%% Activity concentration used to scale sinogram counts with. 
uptakeData   = zeros([vox.petSim.nxn size(data(PTscanNum).data,4)]); %Unattenuated
for j = 1:size(data(PTscanNum).data,4)
    for i = 1:vox.petSim.nxn(3)
        % pad to match PET
        tmp0  = zeros(vox.petSim.nxn(1:2));
        tmp1  = imresize(data(PTscanNum).data(:,:,i,j),vox.pet.xyz(1)/vox.petSim.xyz(1),'cubic');
        
        xA = max(1,round(( vox.petSim.nxn(1) - size(tmp1,2) )/2) + 1);
        xB = xA + size(tmp1,2) - 1;
        yA = max(1,round(( vox.petSim.nxn(1) - size(tmp1,1) )/2) + 1);
        yB = yA + size(tmp1,1) - 1;
        
        tmp0(yA:yB,xA:xB) = tmp1;
        uptakeData(:,:,i,j) = imresize( tmp0,vox.petSim.nxn(1:2), 'cubic' );
    end
end
uptakeData(uptakeData < 0)     = 0; % Alignment verified
FWuptakeData = ones([vox.petSim.rtz(1:3) size(data(PTscanNum).data,4)]);
for j = 1:size(data(PTscanNum).data,4)
    for i = 1:vox.petSim.rtz(3)
        FWuptakeData(:,:,i,j) = FWAC(:,:,i).*radon( uptakeData(:,:,i,j), PHI );
    end
end
activityConc = sum( reshape(FWuptakeData,[prod(vox.petSim.rtz) size(data(PTscanNum).data,4)]),1 );
activityConc = activityConc/max(activityConc); %normalized

return
