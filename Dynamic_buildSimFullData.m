function  [FWPTtrue,FWPTscatter,FWPTrandoms,FWAC,wcc,counts] = Dynamic_buildSimFullData(...
    refPT,FWAC,PSFsim,scatterK,vox,countScale,SF,RF,sensScale)
%%
%******************************************************************************************************
%"Dynamic_buildSimFullData"
%   Builds un-noised PET projection data and attenuation projections
%
% CRS, 08/01/2013
% IH,  04/07/2016
%
% Usage: 
%   [FWPTtrue,FWPTscatter,FWPTrandoms,FWAC,wcc] = Dynamic_buildSimFullData(refPT,muCT,psf,vox,countsTotal,SF,RF)
%       refPT       = reference PET image
%       FWAC        = forward projection of mumap
%       PSFsim      = Postfilter matched to simulation size
%       scatterK    = Scatter kernel matched to simulation size
%       vox         = voxel sizes, dimensions for input, simulation and output images
%       countsScale = mean total counts per active voxel
%       SF          = Scatter fraction
%       RF          = Randoms fraction
%       sensScale   = Sensitivity scale factor for each slice
%
%******************************************************************************************************
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

%% Builds un-noised data 

% rescale Images to match scanner binning
PTtrue = zeros(vox.petSim.nxn);
for i = 1:vox.petSim.nxn(3)
    PTtrue(:,:,i)  = imresize(refPT(:,:,i),vox.petSim.nxn(1:2),'cubic')*sensScale(i);
end
PTtrue(PTtrue < 0) = 0;

% scale PET for counts
cntsPET   = PTtrue; cntsPET(cntsPET <= 0.1) = 0; cntsPET(cntsPET > 0) = 1;
countsPET = countScale * vox.pet.vol/1000 * (vox.pet.nxn(1)/vox.petSim.nxn(1))^2 * ...
    sum(cntsPET(:));

% Blur images
PTscatter = zeros(size(PTtrue));
for i = 1:vox.petSim.nxn(3)
    PTtrue(:,:,i)    = imfilter(PTtrue(:,:,i), PSFsim,      'replicate','same','conv');
    PTscatter(:,:,i) = imfilter(PTtrue(:,:,i), scatterK, 'replicate','same','conv');
end

%% Forward project images
PHI         = 0:180/vox.petSim.rtz(2):180*(1-1/vox.petSim.rtz(2));
FWPTtrueNAC = zeros([vox.petSim.rtz(1) vox.petSim.rtz(2) vox.petSim.rtz(3)]);
FWPTscatter = zeros([vox.petSim.rtz(1) vox.petSim.rtz(2) vox.petSim.rtz(3)]);
FWPTrandoms = zeros([vox.petSim.rtz(1) vox.petSim.rtz(2) vox.petSim.rtz(3)]);
for i = 1:vox.petSim.rtz(3)
    FWPTtrueNAC(:,:,i) = radon( PTtrue(:,:,i),    PHI );
    FWPTscatter(:,:,i) = radon( PTscatter(:,:,i), PHI );
end
radA = floor( (vox.petSim.rtz(1) - vox.petSim.nxn(1))/2 ) + 1; radB = vox.petSim.rtz(1) - radA - 1;
FWPTrandoms(radA:radB,:,:) = 1;
FWPTtrue = FWPTtrueNAC.*FWAC;

% Attenuate total counts
countsTrue    = countsPET * sum(FWPTtrue(:)) / sum(FWPTtrueNAC(:));
if isnan(countsTrue); countsTrue = 0; end %IDA
countsScatter = SF/(1-SF) * countsTrue;
countsRandoms = RF/(1-RF) * (countsTrue + countsScatter);

wcc = countsTrue / sum( FWPTtrue(:) );
wcc = wcc * (  vox.petSim.nxn(1) / vox.petOut.nxn(1) ); 
wcc = wcc*sensScale;
wcc(wcc == 0)  = Inf; %IDA
wcc(isnan(wcc))= Inf; %IDA

FWPTtrue    = countsTrue    * FWPTtrue    / sum( FWPTtrue(:)    );
FWPTscatter = countsScatter * FWPTscatter / sum( FWPTscatter(:) );
FWPTrandoms = countsRandoms * FWPTrandoms / sum( FWPTrandoms(:) );
FWPTtrue(isnan(FWPTtrue))       = 0; %IDA
FWPTscatter(isnan(FWPTscatter)) = 0; %IDA
FWPTrandoms(isnan(FWPTrandoms)) = 0; %IDA

counts.total   = countsTrue + countsScatter + countsRandoms;
counts.true    = countsTrue;
counts.scatter = countsScatter;
counts.randoms = countsRandoms;
counts.NEC     = counts.true^2 / counts.total;
counts.ID      = counts.NEC / prod(vox.petOut.nxn);

counts
end
