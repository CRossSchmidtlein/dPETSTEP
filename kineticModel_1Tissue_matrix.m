function Cpet = kineticModel_1Tissue_matrix(t,dt,Cp,p)
%%
%*******************************************************************************************************
%| Generates a response TAC corresponding to the 1-tissue compartment model from a given input of      |
%| kinetic parameters, frame vector and arterial input function.                                       |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| np = number of model parameters = 2.                                                                |
%| nv = number of voxels.                                                                              |
%|                                                                                                     |
%| USAGE  :  Cpet = kineticModel_1TissueModel_matrix(t,dt,Cp,p).                                       |
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
%|                          E.g. for the 1-tissue model: param(v,:) = [K1 k2].                         |
%|                          p = [p_1,1 p_1,2 ... p_1,np ]                                              |
%|                              [ ...      ...     ...  ]                                              |
%|                              [p_nv,1    ...   p_nv,np].                                             |
%|                                                                                                     |
%| OUTPUT :  Cpet           Matrix of response tissue TAC, size [nv,(f-1)], unit (same as Cp).         |
%|                          Cpet = [Cpet_1,1 Cpet_1,2 ... Cpet_1,f-1 ]                                 |
%|                                 [Cpet_2,1 Cpet_2,2 ... Cpet_2,f-1 ]                                 |
%|                                 [   ...            ... Cpet_nv,f-1].                                |
%|               __________________                                                                    |
%|      |    |  | K1   _________   |                                                                   |
%|      | Cp |--|---->|         |  |                                                                   |
%|      |    |<-|-----|   C1    |  |                                                                   |
%|      |    |  | k2  |_________|  |                                                                   |
%|      |    |  |__________________| Cpet                                                              |
%|                                                                                                     |
%|   Theoretical 1-tissue model:                                                                       |
%|   *-------------------------------------------*                                                     |
%|   | Cpet   = [ K1*exp(-k2*t) ] CONV [Cp]      |                                                     |
%|   *-------------------------------------------*                                                     |
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
term1           = bsxfun( @times, -p(:,2), t );           %-k2*t
C               = bsxfun( @times, p(:,1), exp( term1 ) ); %K1*exp(term1)
% fprintf('Done!\n')

%% Remove NaN and Inf values.
C(isnan(C))     = 0;
C(isinf(C))     = 0;

%% Convolve with blood to get final measured tissue TAC.
cellC  = num2cell(C, 2); % creates a cell array where each element is a row of A
Cpet   = dt*cell2mat( cellfun(@(row)(conv(row, Cp')), cellC, 'UniformOutput',0) )';
Cpet   = Cpet(1:numel(Cp),:);
    