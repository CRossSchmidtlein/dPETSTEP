function w = calculateWeights(varargin)
%%
%******************************************************************************************************
% Calculates frame weights = 1/variance according to user specified model. 
%
% IH, 19/04/2016
%
% f = number of frames
%
% INPUT:    t        - mid frame (time) column vector [f-1,1], unit (sec).
%           C        - tissue TAC [f-1,1].
%           type     - type of weighting, string 'w1', 'w2'...
%           frameVar - frame variance, column vector [f-1,1].
%           halflife - halflife of nuclide in (sec), scalar.
%
%           Available weight model types are:
%           'ones' : uniform weighting (vector of ones)
%           'w1'   : 1/frameVariance
%           'w2'   : decayFactor^2
%           'w3'   : decayFactor * 1/(frameLength*TAC )
%           'w4'   : frameLength * exp(-lambda*t) / TAC
%           'w5'   : frameLength * exp(-lambda*t)
%
% OUTPUT:   w        - weight vector [f-1,1].
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

%% Parse input.
for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 't'
            t           = varargin{i+1};  % Mid frame (time) column vector [f,1].
        case 'C'
            C           = varargin{i+1};  % Tissue TAC [f,1].
        case 'type'
            type        = varargin{i+1};  % Type of weighting, string 'ones','w1', 'w2'...
        case 'frameVar'
            frameVar    = varargin{i+1};  % Frame variance [f,1].
        case 'halflife'
            halflife    = varargin{i+1};  % Nuclide halflife in sec, scalar.
        otherwise
            fprintf('Unknown argument ''%s''.\nExiting...\n',varargin{i});
            w = [];
            return;
    end
end

%% Check input.
if exist('frameVar','var')
    if size(frameVar,2) > size(frameVar,1)
        frameVar = frameVar';
    end
elseif strcmp(type,'w1')
    fprintf('You have to supply the frame variance for weight ''w1''! Exiting...')
    w = [];
    return;
end

%% Column vector.
if size(t,2) > size(t,1)
    t = t';
end
if exist('C','var')
    if size(C,2) > size(C,1)
        C = C';
    end
end

%% Frame length. Assume first frame from t=0.
frameLength     = zeros(size(t));
frameLength(1)  = t(1)*2;
for i = 2:numel( t )
    low              = sum( frameLength(1:i-1) );
    frameLength(i,1) = 2*( t(i) - low );
end
        
%% Decay constant and decay factor for some weight types.
switch type
    case{'w2','w3','w4','w5'}
        if ~exist('halflife','var')
            fprintf('No nuclide halflife specified!\n')
            w = [];
            return;
        end
        % Decay factor
        lambda          = log(2)/halflife; %1/sec
end
switch type
    case {'w2','w3'}        
        % A_decay = decayFactor * A_noDecay
        decayFactor     = zeros( size(t));
        decayFactor(1)  = ( 1 - exp(-lambda*frameLength(1)) )...
            ./ ( lambda*frameLength(1) );
        for i = 2:numel( t )
            decayFactor(i,1) = (exp(-lambda*(t(i)-frameLength(i)/2)) - exp(-lambda*(t(i)+frameLength(i)/2)) )...
                ./ ( lambda*frameLength(i) );
        end
end

%% Weights, ideally proportional to 1/sigma^2
switch type
    case 'ones'
        w           = ones(size(t));
    case 'w1'
        w           = frameVar.^(-1);
    case 'w2'
        w           = decayFactor.^(2);
    case 'w3'
        w           = decayFactor .* (frameLength .* C ).^(-1);
    case 'w4'
        w           = frameLength .* exp(-lambda*t) .* ( C.^(-1) );
    case 'w5'
        w           = frameLength .* exp(-lambda*t);
end

%% Remove inf values and normalize to 1
%     infIndex    = isinf(w);
infIndex    = w > 1e4;
w(infIndex) = min( w( logical(1-infIndex) ) );
w           = w/sum( w );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%