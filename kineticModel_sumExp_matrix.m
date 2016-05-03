function Cpet = kineticModel_sumExp_matrix(t,dt,Cp,p)
%%
%*******************************************************************************************************
%| Generates a response TAC corresponding to an arbitrary sum of exponentials, from a given input of   |
%| kinetic parameters, frame vector and arterial input function.                                       |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| np = number of model parameters.                                                                    |
%| nv = number of voxels.                                                                              |
%| N  = number of exponentials.                                                                        |
%|                                                                                                     |
%| USAGE  :  Cpet = kineticModel_sumExp_matrix(t,dt,Cp,p).                                             |
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
%|                          p = [p_1,1 p_1,2 ... p_1,np ]                                              |
%|                              [ ...      ...     ...  ]                                              |
%|                              [p_nv,1    ...   p_nv,np].                                             |
%|                                                                                                     |
%| OUTPUT :  Cpet           Matrix of response tissue TAC, size [nv,(f-1)], unit (same as Cp).         |
%|                          Cpet = [Cpet_1,1 Cpet_1,2 ... Cpet_1,f-1 ]                                 |
%|                                 [Cpet_2,1 Cpet_2,2 ... Cpet_2,f-1 ]                                 |
%|                                 [   ...            ... Cpet_nv,f-1].                                |
%|      |    |   ___________________                                                                   |
%|      | Cp |  |  p1   _________   |                                                                  |
%|      |    |--|----->|         |  |                                                                  |
%|      |    |<-|------|   C1    |  |                                                                  |
%|      |    |  |  p2  |_________|  |                                                                  |
%|      |    |  |  p3   _________   |                                                                  |
%|      |    |--|----->|         |  |                                                                  |
%|      |    |<-|------|   C2    |  |                                                                  |
%|      |    |  |  p4  |_________|  |                                                                  |
%|      |    |  |  .        .       |                                                                  |
%|      |    |  |  .        .       | Cpet                                                             |
%|      .    .  .  .        .       .                                                                  |
%|                                                                                                     |
%|   Theoretical sum of exponentials model:                                                            |
%|   *----------------------------------------------------------------------------------*              |
%|   | Cpet =  [ p1*exp(-p2*t) + p3*exp(-p4*t) + ... + p2N-1*exp(-p2N*t) ] CONV [ Cp ]  |              |
%|   *----------------------------------------------------------------------------------*              |
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

%% No of parameters (= 2*number of exponentials).
NP = floor(size(p,2)/2)*2; %Even number, exclude spillover even if Vp specified

%%  C = p1*exp(-p2*t) + p3*exp(-p4*t) + ... + pNP-1*exp(-pNP*t).
% fprintf('Calculating sum of exponentials... ')
C = 0;
for i = 1:2:NP
    expterm = exp( bsxfun( @times, -p(:,i+1),t ) );
    C       = bsxfun( @times, p(:,i), expterm ) + C;  
end
clear expterm
% fprintf('Done!\n')

%% Remove NaN and Inf values.
C(isnan(C)) = 0;
C(isinf(C)) = 0;

%% Convolve with blood to get final measured tissue TAC, Cpet = [C] CONV [Cp].
% fprintf('Convolving C with Cp... ')
cellC       = num2cell(C, 2); % creates a cell array where each element is a row of C
Cpet        = dt*cell2mat( cellfun(@(row)(conv(row, Cp')), cellC, 'UniformOutput',0) )';
Cpet        = Cpet(1:numel(Cp),:);
% fprintf('Done!\n')  
