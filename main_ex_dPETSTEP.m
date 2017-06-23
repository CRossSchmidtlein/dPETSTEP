%% Settings
% Settings will be read from "Dynamic_setSimParameters.m"
% Make sure flag "simSet.FBP_OUT" is set to "true".

% Load example data
load data_ex_dPETSTEP.mat

%% Run simulation
[data,simSet,FBP4D,OS4D,OSpsf4D,counts,countsNoise,nFWprompts,FWtrues,FWscatters,FWrandoms,wcc] = Dynamic_main(data,frame,Cif,scaleFactor);

%% Model fitting
p0         = 0.001*ones(1,5);
paramImage = modelFitting_main('image',FBP4D,...                                      
                                   'model','2Tissue',... 
                                   'p0',p0,...
                                   'midFrame',midFrame,...                                                                 
                                   'Cp',Cif);                               