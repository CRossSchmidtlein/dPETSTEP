function mumap = convertCT2mumap(CT)
%%
%******************************************************************************************************
% Converts the HU values of the CT to a mu-map (1/cm).
% Reference: C Burger et al., "PET attenuation coefficients from CT images: 
% Experimental evaluation of the transformation of CT into PET 511-keV 
% attenuation coefficients", EJNM 29(7), pp.922-927 (2002).
%
% IH, 19/04/2016
% 
% USAGE  :  mumap = convertCT2mumap(CT).
%
% INPUT  :  CT      Matrix in Hounsfield units, size [x,y,z].
%
% OUTPUT :  mumap   Matrix in units (1/cm), size [x,y,z].
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

%% Attenuation coeff at different energies, unit 1/cm
mu_water_511keV = 0.096;
mu_bone_511keV  = 0.172;
mu_water_80keV  = 0.184;
mu_bone_80keV   = 0.428;

%% Allocate 
mumap           = zeros(size(CT));

%% Convert
indx            = find(CT <= 0);
mumap(indx)     = mu_water_511keV * ( CT(indx)+1000 )/1000;
indx            = find(CT > 0);
mumap(indx)     = mu_water_511keV + CT(indx)*mu_water_80keV* ...
                  (mu_bone_511keV - mu_water_511keV) / (1000*(mu_bone_80keV - mu_water_80keV));
mumap(mumap<0)  = 0;

end