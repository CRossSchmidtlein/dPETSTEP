function OS = recon_osemPSFSimData_nonUniformSliceSens(nFWPTtotal,RS,CTAC,initPT,PSF, ...
    wcc,vox,iterNUM,subNUM)
%%
%******************************************************************************************************
%"recon_osemPSFSimData_nonUniformSliceSens"
%   Reconstructs PET-like images via OSEM w/ psf for non-uniform axial sensitivity
%
% CRS, 08/01/2013
% IH,  19/04/2016
%
%Usage:
%   [OS] = recon_osemPSFSimData_nonUniformSliceSens(nFWPTtotal, RS,CTAC,PSF,POST,wcc,radBin,tanBin,simSize,zSlice,iterNUM,reconType)
%       nFWPTtotal = Projection data
%       RS         = scatter + random projection reference data
%       CTAC       = CT attenuation correction data
%       PSF        = FWHM of PSF
%       POST       = FWHM of post filter
%       wcc        = "Well-counter correction"
%       radBin     = radial bins
%       tanBin     = projection bins
%       simSize    = Final image size
%       zSlice     = number of slices
%       iterNUM    = number of iterations
%
%******************************************************************************************************
% Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team.
%
% This file is part of The Computational Environment for Radiotherapy Research (CERR).
%
% CERR development has been led by:  Aditya Apte, Divya Khullar, James Alaly, and Joseph O. Deasy.
%
% CERR has been financially supported by the US National Institutes of Health under multiple grants.
%
% CERR is distributed under the terms of the Lesser GNU Public License.
%
%     This version of CERR is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
% CERR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CERR.  If not, see <http://www.gnu.org/licenses/>._nonUniformSliceSens

%% Timing
tic; 

%% Reconstructs OS images w/ psf
supportCoeff = 0.99;
% resize images to odd dimensions
OddMat = true;
if (mod(vox.petOut.nxn(1),2) == 0)
    vox.petOut.nxn(1:2) = vox.petOut.nxn(1:2) + 1;
    
    vox.petOut.rtz(1) = 2*ceil(norm( vox.petOut.nxn(1:2) - floor(( vox.petOut.nxn(1:2)-1 )/2)-1)) + 3;
    
    initPTup = initPT;
    initPT = zeros(vox.petOut.nxn);
    for i = 1:vox.petOut.nxn(3)
        initPT(:,:,i) = imresize(initPTup(:,:,i),[vox.petOut.nxn(1:2)]);
    end
    
    OddMat = false;
end

nn = (vox.petOut.nxn(1)-1)/2;
x2 = linspace(-nn,nn,vox.petOut.nxn(1)).^2;
disk = x2(ones(vox.petOut.nxn(1),1),:) + x2(ones(1,vox.petOut.nxn(1)),:)' <= (vox.petOut.nxn(1)/2)^2;

% extraineous projection padding for corners
radPadA = floor((vox.petOut.rtz(1) - vox.petOut.nxn(1))/2) + 1;
radPadB = radPadA + vox.petOut.nxn(1) - 1;
radRefPadA = floor((vox.petSim.rtz(1) - vox.petSim.nxn(1))/2) + 1;
radRefPadB = radRefPadA + vox.petSim.nxn(1) - 1;

%wcc = wcc*(radRef/rad)^2

% Build subset angles
subTanBin = vox.petOut.rtz(2)/subNUM;
if (mod(subTanBin,1) ~= 0)
    fprintf('ERROR: projection angles not divisible by subsets\n');
    fprintf('\tProjections = %3d, \tSubsets = %2d, \t Number of angle bins = %3.2f\n', ...
        vox.petOut.rtz(2),subNum,subTanBin);
    return;
end

PHIref  = 0:180/vox.petOut.rtz(2):180*(1-1/vox.petOut.rtz(2));
PHI     = zeros([subTanBin subNUM]);
PHIbins = zeros([subTanBin subNUM]);
phiVec  = randperm(subNUM) - 1;
for i = 1:subNUM
    j = phiVec(i);
    PHI(:,i) = j*180/vox.petOut.rtz(2) : subNUM*180/vox.petOut.rtz(2) : ...
        180*(1-(subNUM-j)/vox.petOut.rtz(2));
    PHIbins(:,i) = j+1 : subNUM : vox.petOut.rtz(2)-subNUM+j+1;
end

% Build data/additive noise/normalization subsets
nFW     = zeros([vox.petSim.rtz(1)  subTanBin subNUM vox.petSim.rtz(3)]);
rs      = zeros([vox.petSim.rtz(1)  subTanBin subNUM vox.petSim.rtz(3)]);
subFWAC = zeros([vox.petSim.rtz(1)  subTanBin subNUM vox.petSim.rtz(3)]);
subH1   = zeros([vox.petOut.nxn(1:2) subNUM vox.petOut.nxn(3)]);
for i = 1:vox.petSim.rtz(3)
    for j = 1:subNUM
        nFW(:,:,j,i)     = nFWPTtotal(:,PHIbins(:,j),i);
        rs(:,:,j,i)      = RS(:,PHIbins(:,j),i);
        subFWAC(:,:,j,i) = CTAC(:,PHIbins(:,j),i);
        CTACtmp          = imresize(subFWAC(:,:,j,i),[vox.petOut.rtz(1) subTanBin],'bilinear');
        subHTmp          = iradon(CTACtmp,PHI(:,j),'linear','none',vox.petOut.nxn(1));
        CTACtmp(CTACtmp < 0) = 0;
        CTsuptmp         = CTACtmp;
        CTsuptmp(CTsuptmp >= supportCoeff) = 0; CTsuptmp(CTsuptmp > 0) = 1;
        CTsupport(:,:,j,i) = CTsuptmp;
        subH1(:,:,j,i)   = disk.*imfilter(subHTmp,PSF','replicate','same','conv');        
    end
end
subH1(subH1 < 1E-12) = 1;

clear CTAC nFWPTtotal RS

if (OddMat)
    OS = zeros([vox.petOut.nxn iterNUM]);
else
    OS = zeros([vox.petOut.nxn(1:2)-1 vox.petOut.nxn(3) iterNUM]);
end
hf      = zeros(vox.petSim.rtz(1),subTanBin);
gfRatio = zeros(vox.petOut.rtz(1),subTanBin);
for i = 1:vox.petOut.rtz(3)
%     tic
    f   = initPT(:,:,i)*wcc(i)*subNUM;
    for j = 1:iterNUM
        for k = 1:subNUM
            
            if (numel(PSF) ~= 1)
                
                % Resizes image data to match projection size to perserve
                % Poisson characteristics
                bf = imfilter(f,PSF,'replicate','same','conv');
                hfTmp = radon(bf/subNUM,PHI(:,k));
                hf(radRefPadA:radRefPadB,:) = ...
                    imresize(hfTmp(radPadA:radPadB,:),[vox.petSim.nxn(1) subTanBin],'bilinear') + ...
                    rs(radRefPadA:radRefPadB,:,k,i)./subFWAC(radRefPadA:radRefPadB,:,k,i);
                % Resizes projection data to match image size
                gfRatio(radPadA:radPadB,:) = ...
                    imresize(nFW(radRefPadA:radRefPadB,:,k,i)./hf(radRefPadA:radRefPadB,:), ...
                    [vox.petOut.nxn(1) subTanBin],'bilinear');
                gfRatio = gfRatio.*CTsupport(:,:,k,i);
                hTg = iradon(gfRatio,PHI(:,k),'linear','none',vox.petOut.nxn(1));
                hTg = imfilter(hTg,PSF','replicate','same','conv');
                
            else
                
                % Resizes image data to match projection size to perserve
                % Poisson characteristics
                hfTmp = radon(f/subNUM,PHI(:,k));
                hf(radRefPadA:radRefPadB,:) = ...
                    imresize(hfTmp(radPadA:radPadB,:),[vox.petSim.nxn(1) subTanBin],'bilinear') + ...
                    rs(radRefPadA:radRefPadB,:,k,i)./subFWAC(radRefPadA:radRefPadB,:,k,i);
                % Resizes projection data to match image size
                gfRatio(radPadA:radPadB,:) = ...
                    imresize(nFW(radRefPadA:radRefPadB,:,k,i)./hf(radRefPadA:radRefPadB,:), ...
                    [vox.petOut.nxn(1) subTanBin],'bilinear');
                gfRatio = gfRatio.*CTsupport(:,:,k,i);
                hTg = iradon(gfRatio,PHI(:,k),'linear','none',vox.petOut.nxn(1));
                
            end
            
            f = disk.*f.*hTg./subH1(:,:,k,i);
            
            f(f < 0) = 0;
            f(isnan(f)) = 0;
            f(isinf(f)) = 0;
            
        end
        if (OddMat)
            OS(:,:,i,j) = disk.*f/wcc(i)/subNUM;
        else
            OS(:,:,i,j) = imresize(disk.*f/wcc(i)/subNUM,size(f)-1);
        end

    end
%     toc
end

fprintf('OSEMpsf time: %.2f sec\n',toc);

end