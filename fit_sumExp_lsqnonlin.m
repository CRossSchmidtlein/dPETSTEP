function [p] = fit_sumExp_lsqnonlin(t, C, Cp, w, p0, dt, doPlot, lowerBounds, upperBounds, algorithm)
%%
%*******************************************************************************************************
%| Fits a response TAC to an arbitrary sum of exponentials, from given time sampling, input function   |
%| and frame weights. Generates a vector of fitted parameters.                                         |
%|                                                                                                     |
%| IH, 19/04/2016                                                                                      |
%|                                                                                                     |
%| f  = number of frames.                                                                              |
%| N  = number of exponentials.                                                                        |
%| np = number of model parameters = 2N or 2N+1.                                                       |
%|                                                                                                     |
%| USAGE  :  p = fit_sumExp_lsqnonlin(t,C,Cp,w,p0,dt,doPlot,lowerBounds,upperBounds,algorithm).        |
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
%|                          p0 = [p0_1 p0_2 p0_3...p0_np].                                             |
%|           dt             Scalar with frame duration of t, unit (s).                                 |
%|           doPlot         Flag to plot fitted solution (1) or not (0).                               |
%|           lowerBounds    Vector with lower bounds for estimate of p, size [1,np]. Default zero.     |
%|           upperBounds    Vector with upper bounds for estimate of p, size [1,np]. Default 100*p0.   |
%|           algorithm      String with desired fitting algorithm. Default 'trust-region-reflective'.  |
%|                                                                                                     |
%| OUTPUT :  p              Vector with fitted model parameters, size [1,np].                          |
%|                          p = [p_1 p_2 p_3...p_np].                                                  |
%|      |    __|____________________                                                                   |
%|      |Cp |  |   p1   _________   |                                                                  |
%|      |   |  |------>|         |  |                                                                  |
%|      |   |  |<------|   C1    |  |                                                                  |
%|      |   |  |   p2  |_________|  |                                                                  |
%|      |   |  |   p3   _________   |                                                                  |
%|      |   |  |------>|         |  |                                                                  |
%|      |   |  |<------|   C2    |  |                                                                  |
%|      |   |  |   p4  |_________|  |                                                                  |
%|      | p(2N+1)  .        .       |                                                                  |
%|      |   |  |   .        .       | Cpet                                                             |
%|      .   .  .   .        .       .                                                                  |
%|                                                                                                     |
%|   Theoretical model:                                                                                |
%|   C        = [ p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t) + ... + p(2N-1)*exp(-p(2N)*t) ] CONV [ Cp ]    |
%|   Cpet     = (1-p(2N+1))*C + p(2N+1)*Cp                                                             |
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

%% Number of exponentials.
noExp         = floor(numel(p0)/2);
%If spillover term (last element of p).
if mod(numel(p0),2) ~= 0
    spillover = 1; 
else
    spillover = 0; 
end 

%% Bounds for the model parameters.
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
p             = lsqnonlin(@fitFunction,p0,lowerBounds,upperBounds,options,t,C,Cp,w,dt,spillover,noExp);
 
%% Plot measured and fitted TAC.
if doPlot
    figNo     = 1;
    C_fitted  = myFunction(p,t,Cp,dt,spillover,noExp);
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
function DIFF = fitFunction(p,t,C_meas,Cp,w,dt,spillover,noExp)
% INPUT: p          Parameter vector [p1,p2,p3,...].
%        t          Time vector, unit [sec].
%        C_meas     tissue response TAC.
%        Cp         Input function.
%        spillover  Include spillover (1) or not (0).
%        w          Weights for each frame.
%        noExp      Number of exponentials.

C           = 0;
for i = 1:2:2*noExp
    C       = p(i)*exp(-p(i+1)*t) + C;
end
C           = conv( C, Cp);            % Convolve with Cp
C           = dt.*C( 1:numel(t) );     % Cut away extra elements from convolution

if spillover
    C       = ( 1-p(noExp*2+1) )*C + p(noExp*2+1)*Cp;  % Consider blood in tissue
end

% 2-Norm solver
DIFF        = (C_meas - C).*sqrt(w);   % Weighted function to minimize...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = myFunction(p,t,Cp,dt,spillover,noExp)
% INPUT: p          Parameter vector [p1,p2,p3,...].
%        t          Time vector, unit [sec].
%        C_meas     tissue response TAC.
%        Cp         Input function.
%        spillover  Include spillover (1) or not (0).
%        w          Weights for each frame.
%        noExp      Number of exponentials.

C           = 0;
for i = 1:2:2*noExp
    C       = p(i)*exp(-p(i+1)*t) + C;
end
C           = conv( C, Cp);            % Convolve with Cp
C           = dt.*C( 1:numel(t) );     % Cut away extra elements from convolution

if spillover
    C       = ( 1-p(noExp*2+1) )*C + p(noExp*2+1)*Cp;  % Consider blood in tissue
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%