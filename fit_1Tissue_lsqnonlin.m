function [p,fitInfo] = fit_1Tissue_lsqnonlin(t, C, Cp, w, p0, dt, doPlot, lowerBounds, upperBounds, options)
%%
%*******************************************************************************************************
%| Fits a response TAC to the 1-tissue compartment model, from given time sampling, input function     |
%| and frame weights. Generates a vector of fitted parameters.                                         |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| np = number of model parameters = 2 or 3.                                                           |
%|                                                                                                     |
%| USAGE  :  p = fit_1Tissue_lsqnonlin(t,C,Cp,w,p0,dt,doPlot,lowerBounds,upperBounds,options).         |
%|                                                                                                     |
%| INPUT  :  t              Vector of mid frame times (evenly sampled), size [(f-1),1], unit (s).      |
%|                          t = [t_1; t_2; ... ;t_f-1].                                                |
%|           C              Vector with tissue TAC, size [(f-1),1], unit (arbitrary) e.g. (Bq/cc).     |
%|                          C = [C_1; C_2; ... ;C_f-1].                                                |
%|           Cp             Vector with AIF values corresponding to frame mid times,                   |
%|                          size [(f-1),1], unit (arbitrary), e.g. (Bq/cc).                            |
%|                          Cp = [Cp_1; Cp_2; ... ;Cp_f-1].                                            |
%|           w              Vector with frame weights, size [(f-1),1].                                 |
%|                          w = [w_1; w_2; ... ;w_f-1].                                                |
%|           p0             Vector with initial parameter guesses, unit rate const. (1/s). Size [1,np].|
%|                          p0 = [K1_0 k2_0 (Vp_0)].                                                   |
%|           dt             Scalar with frame duration of t, unit (s).                                 |
%|           doPlot         Flag to plot fitted solution (1) or not (0).                               |
%|           lowerBounds    Vector with lower bounds for estimate of p, size [1,np]. Default zero.     |
%|           upperBounds    Vector with upper bounds for estimate of p, size [1,np]. Default 100*p0.   |
%|           options        Structure with options for fitting (from optimset). Empty-->using defaults.|
%|                          options.MaxIter = 1000, options.Algorithm = 'trust-region-reflective',...  |
%|                                                                                                     |
%| OUTPUT :  p              Vector with fitted model parameters, size [1,np].                          |
%|                          p = [K1 k2 (Vp)].                                                          |
%|           fitInfo        Struct with resnorm,residual,exitflag,output,lambda,jacobian from fit.     |
%|      |    __|____________________                                                                   |
%|      |Cp |  |   K1   _________   |                                                                  |
%|      |   |  |------>|         |  |                                                                  |
%|      |   |  |<------|   C1    |  |                                                                  |
%|      |   |Vp|   k2  |_________|  |                                                                  |
%|      |   |__|____________________| Cpet                                                             |
%|      |      |                                                                                       |
%|                                                                                                     |
%|   Theoretical 1-tissue model:                                                                       |
%|   *----------------------------------------*                                                        |
%|   | C      = [ K1*exp(-k2*t) ] CONV [Cp]   |                                                        |
%|   | Cpet   = (1-Vp)*C + Vp*Cp              |                                                        |
%|   *----------------------------------------*                                                        |
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

%% If spillover term Vp (last element of p).
if mod(numel(p0),2) ~= 0
    spillover = 1; 
else
    spillover = 0; 
end 

%% Options
if isempty(options)
    options = optimset('Display','off');
end

%% Non-linear fitting.
fitInfo = struct;
[p,fitInfo.resnorm,fitInfo.residual,fitInfo.exitflag,fitInfo.output,fitInfo.lambda,fitInfo.jacobian] = ...
    lsqnonlin(@fitFunction,p0,lowerBounds,upperBounds,options,t,C,Cp,w,dt,spillover);
 
%% Plot measured and fitted TAC.
if doPlot
    figNo    = 1;
    C_fitted = myFunction(p,t,Cp,dt,spillover);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DIFF = fitFunction(p,t,C_meas,Cp,w,dt,spillover)
% INPUT: p          Parameter vector [K1, k2, (Vp)], units [1/sec].
%        t          Time vector, unit [sec].
%        C_meas     tissue response TAC.
%        Cp         Input function.
%        w          Weights for each frame.

% Analytical solution
C_part      = p(1)*exp( -p(2)*t );
C_part(isnan(C_part)) = 0;
C           = conv( C_part, Cp);          % Convolve with Cp
C           = dt.*C( 1:numel(t) );        % Cut away extra elements from convolution
if spillover
    C       = ( 1-p(3) )*C + p(3)*Cp;     % Consider blood in tissue
end

% 2-Norm solver
DIFF        = (C_meas - C).*sqrt(w);      % Weighted function to minimize...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = myFunction(p,t,Cp,dt,spillover)
% INPUT: p          Parameter vector [K1, k2, (Vp)], units [1/sec].
%        t          Time vector, unit [sec].
%        C_meas     tissue response TAC.
%        Cp         Input function.
%        w          Weights for each frame.

% Analytical solution
C_part      = p(1)*exp( -p(2)*t );
C_part(isnan(C_part)) = 0;
C           = conv( C_part, Cp);          % Convolve with Cp
C           = dt.*C( 1:numel(t) );        % Cut away extra elements from convolution
if spillover
    C       = ( 1-p(3) )*C + p(3)*Cp;     % Consider blood in tissue
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%