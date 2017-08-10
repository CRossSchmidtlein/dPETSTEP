function gui_saveSettings( handles )

%     funcName = findobj(gcf,'Name',funcName);
    %% Excel filename prompt
    [fileName, pathName,ok] = uiputfile('*.xls', 'Choose a file name'); 
%     [fileName,pathName,ok] = uigetfile('*.xls*');
    
    %% If press cancel
    if ~ok; return; end
    
    %% Print
    fprintf('\nSaving settings in file "%s"...',fileName);
    
    %% Initialize settings
    settings = struct;
    
    %% Extract values from edit boxes
    names = get(findobj(gcf,'style','edit'),'Tag');
    for i = 1:numel(names)
        settings.(names{i}) = handles.(names{i}).String;
    end
    
    %% Extract values from radio buttons   
    names = get(findobj(gcf,'style','radiobutton'),'Tag');
    for i = 1:numel(names)
        pushed = get(handles.(names{i}),'value');
        if pushed
            parent                 = handles.(names{i}).Parent;
            settingName            = parent.Tag;
            settings.(settingName) = handles.(names{i}).Tag;
        end
    end
    
    %% Extract values from toggle buttons   
    names = get(findobj(gcf,'style','togglebutton'),'Tag');
    for i = 1:numel(names)
        pushed = get(handles.(names{i}),'value');
        if pushed
            parent                 = handles.(names{i}).Parent;
            settingName            = parent.Tag;
            settings.(settingName) = handles.(names{i}).Tag;
        end
    end
    
    %% Extract values from checkboxes   
    names = get(findobj(gcf,'style','checkbox'),'Tag');
    for i = 1:numel(names)
        settings.(names{i}) = handles.(names{i}).Value;
    end
    
    %% Write to file
    setNames = fieldnames( settings );
    setVals  = struct2cell(settings);
    xlswrite( fullfile(pathName,fileName), [setNames setVals]);

    %% Print
    fprintf(' Done!\n');
end