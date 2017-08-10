function [Cif,ROI_Cif] = getInputFunction(handles) 

Cif     = [];
ROI_Cif = [];

if handles.inputFuncID.Value==1 %image-derived
    if handles.ifROIFromWS.Value==1 %selected read from WS
        ROI_Cif   = evalin('base',handles.inputFuncROIName.String);
    else
        [~,~,ext] = fileparts(handles.inputFuncROIName.String);
        if strcmp(ext,'.mat') %MAT-file
            ROI_Cif = cell2mat(struct2cell( load(handles.inputFuncROIName.String) ));
        else %Text file
            ROI_Cif = load(handles.inputFuncROIName.String);
        end
        
%         fileID = fopen(handles.inputFuncROIName.String,'r');
%         ROI_Cif   = fscanf(fileID,'%f');
%         fclose(fileID);
    end
else %external
    if handles.ifFromWS.Value==1 %selected read from WS
        Cif   = evalin('base',handles.inputFuncName.String);
    else
        [~,~,ext] = fileparts(handles.inputFuncName.String);
        if strcmp(ext,'.mat') %MAT-file
            Cif = cell2mat(struct2cell( load(handles.inputFuncName.String) ));
        else %Text file
            Cif = load(handles.inputFuncName.String);
        end
        
%         fileID = fopen(handles.inputFuncName.String,'r');
%         Cif   = fscanf(fileID,'%f');
%         fclose(fileID);
    end
end

end