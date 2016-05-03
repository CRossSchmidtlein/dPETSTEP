function dynamicImage = createDynamicPETfromParametricImage_matrix(varargin)
%%
%****************************************************************************************************
%| Creates a 4D dynamic image of a 4D matrix with kinetic parameters (parametric image), using      |
%| a given blood input function or reference TAC, and frame vector. PIM = parametric image.         |
%|                                                                                                  |
%| IH, 19/04/2016                                                                                   |
%|                                                                                                  |
%| N = no of kinetic parameters.                                                                    |
%| f = no of frames.                                                                                |
%|                                                                                                  |
%| USAGE  :   dynamicImage = createDynamicPETfromParametricImage_matrix('paramImage',paramImage,... |
%|                                                               'model','2Tissue',...              |
%|                                                               'Cif',[cif_1 cif_2...cif_f-1],...  |
%|                                                               'CifScaling',2,...                 |
%|                                                               'doDecay',halflife,...             |
%|                                                               'frame',[f_1 f_2...f_f],...        |
%|                                                               'dt',0.5,...                       |
%|                                                               'interpMethod',interpMethod).      |
%|                                                                                                  |
%| INPUT  :   paramImage    4D matrix with parametric image - one parameter set per 3D              |
%|                          voxel, unit (1/s) for rate constants.                                   |
%|                          paramImage = [nx*ny*nz*N].                                              |
%|                          paramImage(x,y,z,:) = [N*1].                                            |
%|            model         String of desired kinetic model, e.g. '1Tissue' or '2Tissue'.           |
%|            frame         Vector with frame start and ends, unit (s).                             |
%|                          frame = [f*1].                                                          |
%|                          [frameStart1; frameStart2=frameEnd1; frameStart3=frameEnd2;...].        |
%|            Cif           Vector with input function to model (arterial or reference tissue),     |
%|                          unit (arbitrary), e.g. (Bq/cc).                                         |
%|                          Cif = [(f-1)*1].                                                        |
%|            dt            Scalar with time step length for calculations, unit (s). Smaller        |
%|                          step means a better calculation of the TACs, but a longer               |
%|                          computation time.                                                       |
%|            OPTIONAL :                                                                            |
%|            CifScaling    Scalar to multiply Cif with.                                            |
%|            doParallell   Optional flag to do parallell computing (1) or not (0)(default).        |
%|            doDecay       Optional string to add physical decay to TACs or not. Specify halflife  |
%|                          in (s) or leave out or specify 'none' to not add decay (default).       |
%|            interpMethod  Optional string with interpolation method (interp1). Default 'linear'.  |
%|                                                                                                  |
%| OUTPUT :   dynamicImage  4D matrix with dynamic image (3D image with a 4th time                  |
%|                          dimension), unit (arbitrary, same as Cif (Cp or Cref).                  |
%|                          dynamicImage = [nx*ny*nz*(f-1)].                                        |
%|                                                                                                  |
%****************************************************************************************************
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

%% Start timing.
mainClock = tic;

%% Check that all arguments have a value and that you have enough arguments.
if mod(numel(varargin),2)
    fprintf('Uneven number (%d) of arguments + values: All arguments do not have a value.\nExiting...\n',numel(varargin));
    dynamicImage = [];
    return;
elseif (numel(varargin)/2) < 5
    fprintf('Too few arguments (%d). Must be at least 5 arguments + values.\nExiting...\n',numel(varargin));
    dynamicImage = [];
    return;
end

%% Read all function arguments
for i=1:2:numel(varargin)-1
    switch varargin{i}
        case 'paramImage'
            paramImage = varargin{i+1}; % 3D cell with one parameter set per voxel.
                                        % paramImage{x,y,z} = [K1 k2 ...], unit (1/s).
        case 'frame'
            frame = varargin{i+1};      % Frame start and ends in unit (s).
                                        % frame = [f*1], f = no of frames.
                                        % [frameStart1; frameStart2=frameEnd1; frameStart3=frameEnd2; ...]. 
        case 'Cif'
            Cif = varargin{i+1};        % Input function to the kinetic model, in unit (Bq/cc). 
                                        % Cif = [(f-1)*1], f = no of frames.
        case 'CifScaling'
            CifScaling = varargin{i+1}; % Scale factor for Cif.
        case 'model'
            modelName = varargin{i+1};  % Kinetic model chosen, e.g. '2-tissue'.
        case 'dt'
            dt = varargin{i+1};         % Time step for calculations, unit (s).
        case 'doParallell'
            doParallell = varargin{i+1}; 
        case 'doDecay'
            doDecay = 1;
            halflife = varargin{i+1}; 
        case 'interpMethod'
            interpMethod = varargin{i+1};
        otherwise
            fprintf('Unknown argument ''%s''.\nExiting...\n',varargin{i});
            dynamicImage = []; 
            return;
    end
end

% Check that all required input arguments are given.
if ~exist('frame','var') || ~exist('dt','var') || ~exist('paramImage','var') || ~exist('modelName','var')
    fprintf('Missing arguments!\nExiting...\n',varargin{i});
    dynamicImage = [];
    return;
end

% Check that needed input arguments are given.
if ~exist('Cif','var')
    fprintf('Missing input ''Cif'' for model ''%s''.\nExiting...\n',modelName);
    dynamicImage = [];
    return;
end

% Default values if not set by user.
if ~exist('doParallell','var')
    disp('Use default : do NOT use parallell computing in Matlab.')
    doParallell = 0; %Do not use Matlabs parallell loops.
end
if ~exist('doDecay','var') || strcmp(halflife,'none')
    disp('Use default : do NOT simulate physical decay.')
    doDecay = 0; %Do not add physical decay to TACs.
end
if ~exist('interpMethod','var')
    disp('Use default : Linear interpolation.')
    interpMethod = 'linear';
end

%% Scaling AIF
if exist('CIFscaling','var')
    fprintf('\nScaling Cif (Cp or Cref) with factor ''%.2f''...',CifScaling)
    Cif = Cif*CifScaling;
    fprintf('Done!\n')
end

%% Number of kinetic parameters in parametric image (assume same in all voxels) and sizes.
noP     = size( paramImage,4 ); %no of parameters.
sizePIM = size(paramImage(:,:,:,1)); %image size, excluding time dimension.
if (numel(sizePIM) <= 2); sizePIM(3) = 1; end %2D image.

%% Check kinetic model and that number of parameters is valid for chosen model. 
switch modelName
    case {'1Tissue','1-Tissue','1tissue','1-tissue'}
        modelParamNo = [2 3]; %2 or 3 parameters (without or with blood spillover term).
    case {'2Tissue','2-Tissue','2tissue','2-tissue'}
        modelParamNo = [4 5]; %4 or 5 parameters (without or with blood spillover term).
    case {'FRTM','frtm'}
        modelParamNo = 4;
    case {'SRTM','srtm'}
        modelParamNo = 3;
    case {'sumExp','sumexp','SUMEXP'}
        modelParamNo = []; %arbitrary number
    otherwise
        fprintf('Unknown model ''%s''.\nExiting...\n',modelName);
        dynamicImage = [];
        return;
end

%% Check that all dimensions are correct according to chosen compartment model.
if ~ismember( noP, modelParamNo ) && ~strcmp(modelName,'sumExp'); 
    fprintf(['The parametric image (PIM) contains [%d] kinetic parameters/voxel. ',...
    'The chosen model ''%s'' should have [%s].\nExiting...\n'],noP,modelName,num2str(modelParamNo));
    dynamicImage = [];
    return; 
end

%% Blood spillover term or not.
if ( (~strcmp(modelName,'SRTM') && ~strcmp(modelName,'FRTM')) || ...
        strcmp(modelName,'sumExp') ) && mod(noP,2)~=0
    useSpillover = 1;
    lastPinModel = noP-1; %Do not include last parameter yet (spillover term Vp).
else
    useSpillover = 0;
    lastPinModel = noP;
end

%% Time properties of input function.
midFrame    = frame(1:end-1) + diff(frame)/2;

%% Decay correction
if doDecay
    frameLength = diff(frame);
    lambda      = log(2)/halflife; %1/sec
    decayFactor =( exp(-lambda*frame(1:end-1)) - exp(-lambda*frame(2:end)) ) ./ (lambda*frameLength);
    fprintf('\nWill multiply TACs with decay factors for nuclide with halflife ''%.3e'' sec.\n\n',halflife);
end

%% Interpolate to even, closely spaced time vector.
fprintf('Interpolating TACs to close and evenly spaced vectors (dt=%.2f s)... ',dt);
t_interp                         = (frame(1):dt:frame(end))'; %interpolated time vector, (s).
Cif_interp                       = interp1( midFrame, Cif, t_interp, interpMethod,'extrap'); %interpolated input function.
Cif_interp(isnan(Cif_interp))    = 0; %set NaN-values to zero.
Cif_interp(Cif_interp<0)         = 0;

fprintf('Done!\n\n')

%% Calculate TACs based on parameter sets.
% Only calculate TACs for unique PIM voxels (parameter sets) to save time and memory.
[paramImageUnique,~,uniqueIndex]    = unique( single( reshape(paramImage,[prod(sizePIM),noP]) ),'rows'); % All unique parameter sets, stacked.
clear paramImage

% Values to divide voxel matrix into smaller matrices to loop over. Too large matrix 
% calculations result in memory errors in Matlab, hence the division. 
N        = size(paramImageUnique,1);                    %Total number of voxels to calculate TACs for.
mSize    = 2e3;                                         %Matrix size for calculations.
divs     = floor(N/mSize);                              %Number of divisions to do matrix calculations on.
last     = N-mSize*divs;                                %Last index, if N/mSize not integer.
ind      = reshape((1:N-last)',mSize,[]);               %Indices for each of the 'divs' loops.
ind      = num2cell(ind,1)';                            %Make indices to cell. Each row is the voxel indices to calculate on in loop.
if last ~= 0; ind{ numel(ind)+1 } = (N-last+1:N)'; end; %Last index row, if not even division.
N_ind    = numel(ind);

% Time and/or memory saved by only doing unique voxels, not all voxels.
fprintf( 'Total no of PIM voxels                        : %d\n', prod(sizePIM) );
fprintf( 'No of unique voxels to calculate TACs for     : %d\n', size(paramImageUnique,1) );
fprintf( 'Matrix size (memory and/or time) save factor  : %.3f%%\n', 100*(1-size(paramImageUnique,1)/prod(sizePIM)) );

% Calculate TACs from kinetic parameters, using the resampled time vector.
fprintf( '\nLoop over %.0f matrices of up to %.0f voxels each...\n',numel(ind),mSize );
               
% Check if you should use parallell computing for the loops.
if doParallell %use parallell loops
   
    tmpCpet    = cell(divs,1);
    T_interp   = repmat(t_interp',[numel(ind{1}),1]);
    
    % Start parallell pool
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool;
    end
    
    % Timing.
    timing(1)  = toc(mainClock);
    
    parfor j = 1:numel(ind)
        % Generate TAC from kinetic parameter set.     
        if j < numel(ind)
            switch modelName
                case {'1Tissue','1-Tissue','1tissue','1-tissue'}
                    Cpet_interp    = kineticModel_1Tissue_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'2Tissue','2-Tissue','2tissue','2-tissue'}
                    Cpet_interp    = kineticModel_2Tissue_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'FRTM','frtm'}
                    Cpet_interp    = kineticModel_FRTM_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'SRTM','srtm'}
                    Cpet_interp    = kineticModel_SRTM_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'sumExp','sumexp','SUMEXP'}
                    Cpet_interp    = kineticModel_sumExp_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
            end
        else %last index may not be full
            switch modelName
                case {'1Tissue','1-Tissue','1tissue','1-tissue'}
                    Cpet_interp    = kineticModel_1Tissue_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'2Tissue','2-Tissue','2tissue','2-tissue'}
                    Cpet_interp    = kineticModel_2Tissue_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'FRTM','frtm'}
                    Cpet_interp    = kineticModel_FRTM_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'SRTM','srtm'}
                    Cpet_interp    = kineticModel_SRTM_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'sumExp','sumexp','SUMEXP'}
                    Cpet_interp    = kineticModel_sumExp_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
            end
        end
        
        % Interpolate back to original time sampling.
        tmpCpet{j}                 = interp1(t_interp,Cpet_interp,midFrame,interpMethod,'extrap');

        % Add decay to TAC, if specified
        if doDecay
            tmpCpet{j}             = bsxfun(@times,tmpCpet{j},decayFactor);
        end
        
        toc(mainClock)
    end
    Cpet = [cell2mat(tmpCpet(1:end-1)') cell2mat(tmpCpet(end)')];
    
else %do not use Matlab parallell loop 
    percentage  = 0; % start value of progress
    fprintf('Loop progress : ');
    lstr        = fprintf('0%%   (0 min 0 sec)');
    Cpet = zeros(numel(midFrame),N);

    T_interp    = repmat(t_interp',[numel(ind{1}),1]);
    
    % Timing.
    timing(1)   = toc(mainClock);
    
    for j = 1:numel(ind)
        % Generate TAC from kinetic parameter set.
        if j < numel(ind)
            switch modelName
                case {'1Tissue','1-Tissue','1tissue','1-tissue'}
                    Cpet_interp    = kineticModel_1Tissue_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'2Tissue','2-Tissue','2tissue','2-tissue'}
                    Cpet_interp    = kineticModel_2Tissue_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'FRTM','frtm'}
                    Cpet_interp    = kineticModel_FRTM_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'SRTM','srtm'}
                    Cpet_interp    = kineticModel_SRTM_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},:) );
                case {'sumExp','sumexp','SUMEXP'}
                    Cpet_interp    = kineticModel_sumExp_matrix( T_interp, dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
            end
        else %last index may not be full
            switch modelName
                case {'1Tissue','1-Tissue','1tissue','1-tissue'}
                    Cpet_interp    = kineticModel_1Tissue_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'2Tissue','2-Tissue','2tissue','2-tissue'}
                    Cpet_interp    = kineticModel_2Tissue_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
                case {'FRTM','frtm'}
                    Cpet_interp    = kineticModel_FRTM_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        repmat(Cif_interp,[numel(ind{end}),1]), paramImageUnique(ind{j},:) );
                case {'SRTM','srtm'}
                    Cpet_interp    = kineticModel_SRTM_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        repmat(Cif_interp,[numel(ind{end}),1]), paramImageUnique(ind{j},:) );
                case {'sumExp','sumexp','SUMEXP'}
                    Cpet_interp    = kineticModel_sumExp_matrix( repmat(t_interp',[numel(ind{end}),1]), dt,...
                        Cif_interp, paramImageUnique(ind{j},1:lastPinModel) );
            end
        end
               
        % Interpolate back to original time sampling.
        Cpet(:,ind{j}')            = interp1(t_interp,Cpet_interp,midFrame,interpMethod,'extrap');
        
        % Add decay to TAC, if specified
        if doDecay
            Cpet(:,ind{j}')     = bsxfun(@times, Cpet(:,ind{j}'),decayFactor);
        end
        
        % Display percentage progress.
        if round( j/N_ind*100 ) > percentage
            percentage_new = round( j/N_ind*100 );
            fprintf( repmat('\b',[1,lstr] ) );
            lstr           = fprintf( ['%.0f%%' repmat(' ',[1,4-length(num2str(percentage_new))]) '(%.0f min %.0f sec)'],...
                percentage_new,floor(toc(mainClock)/60),rem(toc(mainClock),60));
            percentage     = percentage_new;
        end
    end
end

%% Remove negative (unphysical) values.
Cpet(Cpet<0) = 0;
        
%% Timing.
timing(2) = toc(mainClock);

%% Clear variables.
clear T_interp Cpet_interp

%% Blood spillover term.
if useSpillover
    fprintf( '\nAdding blood spillover... ' );
    Vp     = paramImageUnique(:,end)';
    Cpet   = bsxfun(@times,(1-Vp),Cpet) + bsxfun(@times,Vp,repmat(Cif,[1,size(Cpet,2)]));
    clear Vp
end
clear paramImageUnique
fprintf('Done!\n')

%% Plug calculated TACs into correct voxels of dynamic image. 
% Reshape dynamic image to same 3D size as parametric image (PIM), but also with 4th time dimension.
fprintf('\nReshaping dynamic image according to PIM... ')
dynamicImage = reshape( Cpet(:,uniqueIndex)',[sizePIM numel(frame)-1] );
fprintf('Done!\n')
clear Cpet

%% End timing.
fprintf('\n--------------------------------------------------\n')
timing(3) = toc(mainClock);
t1        = timing(1)+timing(3)-timing(2); %overhead
t2        = timing(2)-timing(1); %TACs
if ((t1+t2)/60>=1)
    fprintf('Total calculation time : %.3f sec (%.1f min)\n',t2+t1,(t2+t1)/60);
else
    fprintf('Total calculation time : %.3f sec\n',t2+t1);
end
if (t1/60>=1)
    fprintf('	Overhead           : %.3f sec (%.1f min)\n',t1,t1/60);
else
    fprintf('	Overhead           : %.3f sec\n',t1);
end
if (t2/60>=1)
    fprintf('	Dynamic image      : %.3f sec (%.1f min)\n',t2,t2/60);
else
    fprintf('	Dynamic image      : %.3f sec\n',t2);
end
fprintf('	Time/voxel         : %.2e sec (ex. overhead)\n\n',t2/N);

end