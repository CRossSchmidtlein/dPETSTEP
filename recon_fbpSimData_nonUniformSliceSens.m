function FBP = recon_fbpSimData_nonUniformSliceSens(nFWPTtotal,RS,CTAC,wcc,vox)
%%
%******************************************************************************************************
%"recon_fbpSimData_nonUniformSliceSens"
%   Reconstructs PET-like images via FBP for non-uniform axial sensitivity
%
% CRS, 08/01/2013
% IH,  19/04/2016
%
%Usage: 
%   [FBP] = recon_fbpSimData_nonUniformSliceSens(nFWPTtotal,RS,CTAC,POST,wcc,radBin,tanBin,simSize,zSlice)
%       nFWPTtotal = Projection data
%       RS         = scatter + random projection reference data
%       CTAC       = CT attenuation correction data
%       POST       = FWHM of post filter
%       wcc        = "Well-counter correction"
%       vox        = radial bins, tanBin, simSize, zSlice
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
% along with CERR.  If not, see <http://www.gnu.org/licenses/>.

%% Timing
tic; 

%% Reconstructs FBP images
supportCoeff = 0.99;

% resize images to odd dimensions
OddMat = true;
if (mod(vox.petOut.nxn(1),2) == 0)
    vox.petOut.nxn(1:2) = vox.petOut.nxn(1:2) + 1;
    
    vox.petOut.rtz(1) = 2*ceil(norm( vox.petOut.nxn(1:2) - floor(( vox.petOut.nxn(1:2)-1 )/2)-1)) + 3;
    
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

PHI = 0:180/vox.petOut.rtz(2):180*(1-1/vox.petOut.rtz(2));

% Build data/additive noise/normalization subsets
f    = zeros(vox.petOut.nxn);
if (OddMat)
    FBP = zeros(vox.petOut.nxn);
else
    FBP = zeros([vox.petOut.nxn(1:2)-1 vox.petOut.nxn(3)]);
end
CTsupport = CTAC;
CTsupport(CTsupport >= supportCoeff) = 0; CTsupport(CTsupport > 0) = 1;
Data = CTsupport.*(nFWPTtotal - RS)./CTAC;
Data(Data <= 0) = 1E-10;
clear nFW rs FWAC

for i = 1:vox.petOut.rtz(3)
    
    dataTmp = zeros(vox.petOut.rtz(1:2));
    dataTmp(radPadA:radPadB,:) = ...
        imresize(Data(radRefPadA:radRefPadB,:,i), ...
        [vox.petOut.nxn(1) vox.petOut.rtz(2)],'bilinear');
    
    f(:,:,i) = iradon(dataTmp,PHI,'linear','Shepp-Logan',vox.petOut.nxn(1));
    
    if (OddMat)
        FBP(:,:,i) = disk.*f(:,:,i)/wcc(i);
    else
        FBP(:,:,i) = imresize(disk.*f(:,:,i)/wcc(i),size(f(:,:,i))-1);
    end
    
end

fprintf('FBP time: %.2f sec\n',toc);

end
