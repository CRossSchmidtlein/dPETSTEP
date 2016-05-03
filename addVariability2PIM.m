function pim_var = addVariability2PIM(varargin)
%%
%******************************************************************************************************
%| Adds variability to an existing parametric image (PIM). Represents biological diversity.           |
%|                                                                                                    |
%| IH, 19/04/2016                                                                                     |
%|                                                                                                    |
%| m        = number of kinetic parameters.                                                           |
%| nx,ny,nz = image dimensions.                                                                       |
%|                                                                                                    |
%| USAGE  :  pim_var = addVariability2PIM('pim',pim,'scale',10).                                      |
%|                                                                                                    |
%| INPUT  :  pim        2 or 2D parametric image, one vector per voxel. Rate constants in unit        |
%|                      (s^-1).                                                                       |
%|                      pim = [nx*ny*nz*m].                                                           |
%|           scale      Scalar with amount of variability, sigma = pim/scale.                         |
%|                      (higher scale-->less variability).                                            |
%|                                                                                                    |
%| OUTPUT :  pim_var    2 or 2D variablility (~noise) parametric image, one vector per voxel. Rate    |
%|                      constants in unit (s^-1).                                                     |
%|                      pim_var = [nx*ny*nz*m].                                                       |
%|                                                                                                    |
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

%% Read all function arguments.
for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 'pim'
            pim = varargin{i+1};        % 2 or 3D parametric image.
                                        % pim(x,y,z,:) = [K1 k2 ...], unit (1/s).
        case 'scale'
            scale = varargin{i+1};      % Scaling of variability, e.g. 10.
        otherwise
            fprintf('Unknown argument ''%s''.\nExiting...\n',varargin{i});
            pim_var = []; 
            return;
    end
end

%% Number of kinetic parameters.
N = size( pim, numel(size(pim)) );

%% 2D or 3D image. 
% 2 or 3 image dimensions plus one parameter dimension.
pim_var = zeros(size(pim));
if numel(size(pim)) < 4 %2D image, 3rd dim = parameter
    % Add noise to each of the N PIMs.
    for i = 1:N
        sigma            = pim(:,:,i)/scale;
        pim_var(:,:,i)   = normrnd(pim(:,:,:,i),sigma);
    end
else %3D image, 4th dim = parameter
    % Add noise to each of the N PIMs.
    for i = 1:N
        sigma            = pim(:,:,:,i)/scale;
        pim_var(:,:,:,i) = normrnd(pim(:,:,:,i),sigma);
    end
end

% Remove negative values (unphysical).
pim_var(pim_var<0) = 0;

end