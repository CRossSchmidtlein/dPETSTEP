function Cpet = kineticModel_FRTM_matrix(t,dt,Cref,p)
%%
%*******************************************************************************************************
%| Generates a response TAC corresponding to the full reference tissue compartment model (FRTM)        |
%| from a given vector of kinetic parameters, frame vector and reference TAC.                          |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| n  = number of compartments in tissue.                                                              |
%| np = number of model parameters = 4.                                                                |
%| nv = number of voxels.                                                                              |
%|                                                                                                     |
%| USAGE  :  Cpet = kineticModel_FRTM_matrix(t,dt,Cref,p).                                             |
%|                                                                                                     |
%| INPUT  :  t              Matrix of mid frame times (evenly sampled), size [nv,(f-1)], unit (s).     |
%|                          t = [t_1,1 t_1,2 ... t_1,f-1 ]                                             |
%|                              [ ...      ...     ...   ]                                             |
%|                              [t_nv,1    ...   t_nv,f-1].                                            |
%|           dt             Scalar with frame duration of t, unit (s).                                 |
%|           Cref           Vector with reference tissue TAC corresponding to frame mid times,         |
%|                          size [(f-1),1], unit (arbitrary), e.g. (Bq/cc).                            |
%|                          Cref = [Cref_1;Cref_2; ... ;Cref_f-1].                                     |
%|           p              Matrix of model parameter values, size [nv,np], unit rate constants (1/s). |
%|                          E.g. for the FRTM: param(v,:) = [R1 k2 k3 BPnd].                           |
%|                          p = [p_1,1 p_1,2 ... p_1,np ]                                              |
%|                              [ ...      ...     ...  ]                                              |
%|                              [p_nv,1    ...   p_nv,np].                                             |
%|                                                                                                     |
%| OUTPUT :  Cpet           Matrix of response tissue TAC, size [nv,(f-1)], unit (same as Cref).       |
%|                          Cpet = [Cpet_1,1 Cpet_1,2 ... Cpet_1,f-1 ]                                 |
%|                                 [Cpet_2,1 Cpet_2,2 ... Cpet_2,f-1 ]                                 |
%|                                 [   ...            ... Cpet_nv,f-1].                                |
%|    ______________________________________________                                                   |
%|   |                                              |                                                  |
%|   |  |    |   K1   _________   k3    _________   |                                                  |
%|   |  | Cp |------>|         |------>|         |  |                                                  |
%|   |  |    |<------|   C1    |<------|   C2    |  |                                                  |
%|   |  |    |   k2  |_________|  k4   |_________|  |                                                  |
%|   |__|____|______________________________________| Tissue with specific binding, Cpet               |
%|    __|____|____________________                                                                     |
%|   |  |    |   K1'  _________   |                                                                    |
%|   |  |    |------>|         |  |                                                                    |
%|   |  |    |<------|   C1'   |  |                                                                    |
%|   |  |    |   k2' |_________|  |                                                                    |                 
%|   |____________________________| Reference tissue with non-specific binding, Cref                   |
%|                                                                                                     |
%|   Theoretical FRTM:                                                                                 |
%|   *-------------------------------------------------------------------------------*                 |
%|   | R1     = K1/K1'                                                               |                 |
%|   | BPnd   = k3/k4                                                                |                 |
%|   | a      = (k3+k4-c)(c-r)/u                                                     |                 |
%|   | b      = (d-k3-k4)(d-r)/u                                                     |                 |
%|   | c      = (s+u)/2                                                              |                 |
%|   | d      = (s-u)/2                                                              |                 |
%|   | u      = sqrt(s^2-q)                                                          |                 |
%|   | q      = 4k2k4                                                                |                 |
%|   | r      = k2/R1                                                                |                 |
%|   | s      = k2+k3+k4                                                             |                 |
%|   | Cpet   = R1*Cref + R1*( [a*Cref] CONV [exp(-ct)] + [b*Cref] CONV [exp(-dt)] ) |                 |
%|   *-------------------------------------------------------------------------------*                 |
%|   Reference: A.A. Lammertsma et al., Comparison of methods for analysis of clinical                 |
%|   raclopride[11C]studies, J. Cereb. Blood Flow Metab. 16(1), 42-52 (1996).                          |
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
% fprintf('Cplculating individual terms in expression for Cpet... ')
u                   = sqrt( (p(:,2) + p(:,3) + p(:,3)./p(:,4)).^2 - 4*p(:,2).*p(:,3)./p(:,4) );
alpha1              = 0.5*(p(:,2) + p(:,3) + p(:,3)./p(:,4) - u); 
alpha2              = 0.5*(p(:,2) + p(:,3) + p(:,3)./p(:,4) + u); 
a                   = ( p(:,3)+p(:,3)./p(:,4)-alpha1 ).*( alpha1-p(:,2)./p(:,1) )./u;
b                   = ( alpha2-p(:,3)-p(:,3)./p(:,4) ).*( alpha2-p(:,2)./p(:,1) )./u;
term1               = bsxfun( @times, p(:,1), Cref );        %R1*Cref
term2               = bsxfun( @times, p(:,1).*a, Cref );     %R1*a*Cref
term3               = bsxfun( @times, p(:,1).*b, Cref );     %R1*b*Cref
term4               = exp( bsxfun( @times, -alpha1, t) );    %exp( -alpha1*t )
term5               = exp( bsxfun( @times, -alpha2, t) );    %exp( -alpha2*t )
% fprintf('Done!\n')

%% Remove zero, NaN and Inf values.
term1(term1<0) = 0; term1(isnan(term1)) = 0; term1(isinf(term1)) = 0;
term2(term2<0) = 0; term2(isnan(term2)) = 0; term2(isinf(term2)) = 0;
term3(term3<0) = 0; term3(isnan(term3)) = 0; term3(isinf(term3)) = 0;
term4(term4<0) = 0; term4(isnan(term4)) = 0; term4(isinf(term4)) = 0;
term5(term5<0) = 0; term5(isnan(term5)) = 0; term5(isinf(term5)) = 0;

%% Convoluted terms.
for i = 1:size(term2,2)
    convTmp24(:,i)  = conv(term2(:,i),term4(:,i));
    convTmp35(:,i)  = conv(term3(:,i),term5(:,i));
end
convTerm24          = dt*convTmp24(1:numel(Cref),:);
convTerm35          = dt*convTmp24(1:numel(Cref),:);

%% Final measured tissue TAC.
% fprintf('Cplculating theoretical Cpet... ')
Cpet                = term1 + convTerm24 + convTerm35;
Cpet(Cpet<0)        = 0;
% fprintf('Done!\n')
    