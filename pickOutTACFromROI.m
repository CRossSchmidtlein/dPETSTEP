function [TACmatrix,TACvariance] = pickOutTACFromROI( image, roi, type )
%%
%******************************************************************************************************
% IH, 19/04/2016.
%
% USAGE :   [TACmatrix,TACvariance] = pickOutTACFromROI( image, roi, type ).
% 
% INPUT :   image        Matrix with 4D image, [nx,ny,nz,f].
%           roi          Matrix with 3D roi, [nx,ny,nz].
%           type         string with desired type, "voxel" or "average".
%                           average - average TAC in ROI, [f,1].
%                           voxel   - each voxel TAC in ROI, [f,noVox].
% 
% OUTPUT:   TACmatrix    Matrix with TAC from the ROI, average or each voxel, [f,1] or [f,noVox].
%           TACvariance  Matrix with ROI variance, [f,1].
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

% 3D index of ROI
%---------------------------------------
ROIindex    = find( roi(:) > 0 );

% Check if you should pick out TAC for each voxel in ROI, or average in ROI
switch type
    case 'average'
        % Allocate TAC matrix
        %---------------------------------------
        TACmatrix   = zeros( size(image,4), 1 );
        TACvariance = zeros( size(image,4), 1 );
        
        % Loop over frames
        %---------------------------------------
        for i = 1:size(image,4)
            temp             = image(:,:,:,i);
            TACmatrix(i,:)   = mean( temp(ROIindex) );
            TACvariance(i,:) = var( temp(ROIindex) );
        end
        
    case 'voxel'
        % Allocate TAC matrix
        %---------------------------------------
        TACmatrix   = zeros( size(image,4), numel(ROIindex) );
        TACvariance = zeros( size(image,4), 1 );
        
        % Loop over frames
        %---------------------------------------
        for i = 1:size(image,4)
            temp             = image(:,:,:,i);
            TACmatrix(i,:)   = temp(ROIindex);
            TACvariance(i,:) = var( temp(ROIindex) );
        end

    otherwise
        disp('Faulty type input! Should be "average" or "voxel".');
        TACmatrix=[];
end

end