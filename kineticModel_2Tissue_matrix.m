function Cpet = kineticModel_2Tissue_matrix(t,dt,Cp,p)
%%
%*******************************************************************************************************
%| Generates a response TAC corresponding to the 2-tissue compartment model from a given input of      |
%| kinetic parameters, frame vector and arterial input function.                                       |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| np = number of model parameters = 4.                                                                |
%| nv = number of voxels.                                                                              |
%|                                                                                                     |
%| USAGE  :  Cpet = kineticModel_2TissueModel_matrix(t,dt,Cp,p).                                       |
%|                                                                                                     |
%| INPUT  :  t              Matrix of mid frame times (evenly sampled), size [nv,(f-1)], unit (s).     |
%|                          t = [t_1,1 t_1,2 ... t_1,f-1 ]                                             |
%|                              [ ...      ...     ...   ]                                             |
%|                              [t_nv,1    ...   t_nv,f-1].                                            |
%|           dt             Scalar with frame duration of t, unit (s).                                 |
%|           Cp             Vector with AIF values corresponding to frame mid times,                   |
%|                          size [(f-1),1], unit (arbitrary), e.g. (Bq/cc).                            |
%|                          Cp = [Cp_1;Cp_2; ... ;Cp_f-1].                                             |
%|           p              Matrix of model parameter values, size [nv,np], unit rate constants (1/s). |
%|                          E.g. for the 2-tissue model: p(v,:) = [K1 k2 k3 k4].                       |
%|                          p = [p_1,1 p_1,2 ... p_1,np ]                                              |
%|                              [ ...      ...     ...  ]                                              |
%|                              [p_nv,1    ...   p_nv,np].                                             |
%|                                                                                                     |
%| OUTPUT :  Cpet           Matrix of response tissue TAC, size [nv,(f-1)], unit (same as Cp).         |
%|                          Cpet = [Cpet_1,1 Cpet_1,2 ... Cpet_1,f-1 ]                                 |
%|                                 [Cpet_2,1 Cpet_2,2 ... Cpet_2,f-1 ]                                 |
%|                                 [   ...            ... Cpet_nv,f-1].                                |
%|               ____________________________________                                                  |
%|      |    |  | K1   _________   k3    _________   |                                                 |
%|      | Cp |--|---->|         |------>|         |  |                                                 |
%|      |    |<-|-----|   C1    |<------|   C2    |  |                                                 |
%|      |    |  | k2  |_________|  k4   |_________|  |                                                 |
%|      |    |  |____________________________________| Cpet                                            |
%|                                                                                                     |
%|   Theoretical 2-tissue model:                                                                       |
%|   *------------------------------------------------------------------*                              |
%|   | alpha1 = 0.5*( k2+k3+k4 - sqrt( [k2+k3+k4]^2 - 4*k2*k4) )        |                              |
%|   | alpha2 = 0.5*( k2+k3+k4 + sqrt( [k2+k3+k4]^2 - 4*k2*k4) )        |                              |
%|   | Cpet   = K1/(alpha2-alpha1)*[ (k3+k4-alpha1)*exp(-alpha1*t) +    |                              |
%|   |                (alpha2-k3-k4)*exp(-alpha2*t) ] CONV [ Cp ]       |                              |
%|   *------------------------------------------------------------------*                              |
%|   Note! Blood spillover term is excluded here.                                                      |
%|                                                                                                     |
%*******************************************************************************************************
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

%% Theoretical solution.
% fprintf('Calculating individual terms in expression for Cpet... ')
alpha(:,1)      = 0.5*( p(:,2)+p(:,3)+p(:,4) - sqrt( (p(:,2)+p(:,3)+p(:,4)).^2 - 4*p(:,2).*p(:,4) ) ); %alpha1
alpha(:,2)      = 0.5*( p(:,2)+p(:,3)+p(:,4) + sqrt( (p(:,2)+p(:,3)+p(:,4)).^2 - 4*p(:,2).*p(:,4) ) ); %alpha2
term1           = bsxfun( @times, -alpha(:,1),t );                                                                                     %-alpha1*t
e1              = bsxfun( @times, p(:,1)./( alpha(:,2)-alpha(:,1) ) .* ( p(:,3)+p(:,4)-alpha(:,1) ), exp(term1) );          %K1/(alpha2-alpha1)*(k3+k4-alpha1)*exp(term1)
clear term1
term2           = bsxfun( @times, -alpha(:,2),t );                                                                                     %%-alpha2*t
e2              = bsxfun( @times, p(:,1)./( alpha(:,2)-alpha(:,1) ) .* ( alpha(:,2)-p(:,3)-p(:,4) ), exp(term2) );          %K1/(alpha2-alpha1)*(alpha2-k3-k4)*exp(term2)
clear term2
C                = e1 + e2;
% fprintf('Done!\n')
clear e1 e2

%% Remove NaN and Inf values.
C(isnan(C))     = 0;
C(isinf(C))     = 0;

%% Convolve with blood to get final measured tissue TAC.
% fprintf('Convolving C with Cp... ')
cellC  = num2cell(C, 2); % creates a cell array where each element is a row of C
Cpet   = dt*cell2mat( cellfun(@(row)(conv(row, Cp')), cellC, 'UniformOutput',0) )';
Cpet   = Cpet(1:numel(Cp),:);
% fprintf('Done!\n')

