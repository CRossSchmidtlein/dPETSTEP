function [p] = fit_FRTM_lsqnonlin(t, C, Cref, w, p0, dt, doPlot, lowerBounds, upperBounds, algorithm)
%%
%*******************************************************************************************************
%| Fits a response TAC to the FRTM compartment model, from given time sampling, reference tissue TAC   |
%| and frame weights. Generates a vector of fitted parameters.                                         |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| np = number of model parameters = 4.                                                                |
%|                                                                                                     |
%| USAGE  :  p = fit_FRTM_lsqnonlin(t,C,Cref,w,p0,dt,doPlot,lowerBounds,upperBounds,algorithm).        |
%|                                                                                                     |
%| INPUT  :  t              Vector of mid frame times (evenly sampled), size [(f-1),1], unit (s).      |
%|                          t = [t_1; t_2; ... ;t_f-1].                                                |
%|           C              Vector with tissue TAC, size [(f-1),1], unit (arbitrary) e.g. (Bq/cc).     |
%|                          C = [C_1; C_2; ... ;C_f-1].                                                |
%|           Cref           Vector with reference tissue TAC corresponding to frame mid times,         |
%|                          size [(f-1),1], unit (arbitrary), e.g. (Bq/cc).                            |
%|                          Cref = [Cref_1; Cref_2; ... ;Cref_f-1].                                    |
%|           w              Vector with frame weights, size [(f-1),1].                                 |
%|                          w = [w_1; w_2; ... ;w_f-1].                                                |
%|           p0             Vector with initial parameter guesses, unit rate const. (1/s). Size [1,np].|
%|                          p0 = [R1_0 k2_0 k3_0 BPnd_0].                                              |
%|           dt             Scalar with frame duration of t, unit (s).                                 |
%|           doPlot         Flag to plot fitted solution (1) or not (0).                               |
%|           lowerBounds    Vector with lower bounds for estimate of p, size [1,np]. Default zero.     |
%|           upperBounds    Vector with upper bounds for estimate of p, size [1,np]. Default 100*p0.   |
%|           algorithm      String with desired fitting algorithm. Default 'trust-region-reflective'.  |
%|                                                                                                     |
%| OUTPUT :  p              Vector with fitted model parameters, size [1,np].                          |
%|                          p = [R1 k2 k3 BPnd].                                                       |
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

%% Bounds for the model parameters [R1 k2 k3 BPnd].
if isempty(lowerBounds)
    lowerBounds = zeros(size(p0));
end
if isempty(upperBounds)
    upperBounds = p0*100;
end

%% Algorithm for fit
if isempty(algorithm)
    algorithm = 'trust-region-reflective';
end

%% Desired fitting options.
options     = optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','off','TolFun',1e-8,'Algorithm',algorithm);

%% Non-linear fitting.
p = lsqnonlin(@fitFunction,p0,lowerBounds,upperBounds,options,t,C,Cref,w,dt);
 
%% Plot measured and fitted TAC.
if doPlot
    figNo    = 1;
    C_fitted = myFunction(p,t,Cref,dt);
    if isempty(figNo)
        figure;
    else
        figure(figNo);
    end
    plot(t/60,C/1000,'o','MarkerSize', 4, 'MarkerFaceColor','b', 'MarkerEdgeColor', 'none'); hold on;
    plot(t/60,C_fitted/1000,'--r','LineWidth', 2); 
    legend('Measured','Best fit','Location','Best')
    hold off;
    xlabel('Time (min)', 'FontSize', 12); ylabel('Activity (kBq*cc^{-1})', 'FontSize', 12);
    xlim([t(1)-0.05*t(1) t(end)*1.01/60]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DIFF = fitFunction(p,t,C_meas,Cref,w,dt)
% INPUT: p          Parameter vector [R1, k2, k3, BPnd], units [1/sec].
%        t          Time vector, unit [sec].
%        C_meas     Tissue response TAC.
%        Cref       Reference tissue TAC.
%        w          Weights for each frame.

% Analytical solution
u            = sqrt( (p(2) + p(3) + p(3)/p(4))^2 - 4*p(2)*p(3)/p(4) );
alpha1       = 0.5*(p(2) + p(3) + p(3)/p(4) - u); 
alpha2       = 0.5*(p(2) + p(3) + p(3)/p(4) + u); 
a            = ( p(3)+p(3)/p(4)-alpha1 )*( alpha1-p(2)/p(1) )/u;
b            = ( alpha2-p(3)-p(3)/p(4) )*( alpha2-p(2)/p(1) )/u;
C1           = dt*conv( a*Cref, exp(-alpha1*t) );  % Convolve with Cref
C2           = dt*conv( b*Cref, exp(-alpha2*t) );  % Convolve with Cref
C            = p(1) + p(1)*( C1(1:numel(t)) + C2(1:numel(t)) ); % Cut away extra elements from convolution

% 2-Norm solver
DIFF         = (C_meas - C).*sqrt(w);      % Weighted function to minimize...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = myFunction(p,t,Cref,dt)
% INPUT: p          Parameter vector [R1, k2, k3, BPnd], units [1/sec].
%        t          Time vector, unit [sec].
%        C_meas     Tissue response TAC.
%        Cref       Reference tissue TAC.
%        w          Weights for each frame.

% Analytical solution
u            = sqrt( (p(2) + p(3) + p(3)/p(4))^2 - 4*p(2)*p(3)/p(4) );
alpha1       = 0.5*(p(2) + p(3) + p(3)/p(4) - u); 
alpha2       = 0.5*(p(2) + p(3) + p(3)/p(4) + u); 
a            = ( p(3)+p(3)/p(4)-alpha1 )*( alpha1-p(2)/p(1) )/u;
b            = ( alpha2-p(3)-p(3)/p(4) )*( alpha2-p(2)/p(1) )/u;
C1           = dt*conv( a*Cref, exp(-alpha1*t) );  % Convolve with Cref
C2           = dt*conv( b*Cref, exp(-alpha2*t) );  % Convolve with Cref
C            = p(1) + p(1)*( C1(1:numel(t)) + C2(1:numel(t)) ); % Cut away extra elements from convolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%