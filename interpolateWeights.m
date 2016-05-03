function w2 = interpolateWeights(t,w,t2)
%******************************************************************************************************
% IH, 19/04/2016
%
% USAGE:    w2 = interpolateWeights(t,w,t2)
%
% INPUT:    t    - mid frame (time) column vector
%           w    - weight vector
%           t2   - mid frame (time) column vector of interpolation time
% 
% OUTPUT:   w2   - interpolated weight vector
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

%% Allocate zero weights
w2            = zeros(size(t2));

%% Divide interpolated weight points by number of points
i  = 1; %original index counter
ii = 1; %new index counter
while ii<=numel(t2)-1 
    if t2(ii)<t(i) && t2(ii+1)>t(i)
        f1 = 1 - abs(t2(ii)-t(i))/abs(t2(ii+1)-t2(ii));
        f2 = 1 - f1;
        w2(ii)   = w(i)*f1;
        w2(ii+1) = w(i)*f2;
        i  = i+1;
        ii = ii+2;
    elseif t2(ii)==t(i)
        w2(ii) = w(i);
        i  = i+1;
        ii = ii+1;
    elseif ii+1 == numel(t2) 
        w2(ii+1) = w(i);
        i  = i+1;
        ii = ii+1;
    else
        ii = ii+1;
    end
end
% figure(1), plot(t,w,'*b',t2,w2,'or')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%