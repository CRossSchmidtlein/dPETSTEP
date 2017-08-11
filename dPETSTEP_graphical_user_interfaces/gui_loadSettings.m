function handles = gui_loadSettings( handles,type )

    %% File calling this file
    switch type
        case 'sim'
            func     = @dPETSTEPgui_sim;
        case 'fit'
            func     = @dPETSTEPgui_fit;
    end
    
    %% Excel filename prompt
    [fileName,pathName,ok] = uigetfile('*.xls*', 'Choose a file');
    
    %% If press cancel
    if ~ok; return; end
    
    %% Print
    fprintf('\nLoading settings from file "%s"...',fileName);
    
    %% Read file into cell
    [~,~,SET]  = xlsread( fullfile(pathName,fileName));
    nanIdx = find( cellfun(@(V) any(isnan(V(:))), SET(:,2)) );
    for i=1:numel(nanIdx); SET{nanIdx(i),2} = ''; end
    
    %% Set values in checkboxes   
    names = get(findobj(gcf,'style','checkbox'),'Tag');
    for i = 1:numel(names)
        indx = find( strcmp( SET(:,1), names{i} ));     
        if ~isempty(indx)
            handles.(names{i}).Value = cell2mat(SET(indx,2));
            func( [handles.(names{i}).Tag '_Callback'], handles.(names{i}),'load',handles)
        end
    end
        
    %% Set values in edit boxes
    names = get(findobj(gcf,'style','edit'),'Tag');
    for i = 1:numel(names)
        indx = find( strcmp( SET(:,1), names{i} ));     
        if ~isempty(indx)
            handles.(names{i}).String = SET{indx,2};
            if ~isempty(SET{indx,2});
                func( [handles.(names{i}).Tag '_Callback'], handles.(names{i}),'load',handles)
            end
        end
    end
    
    %% Set values in radio buttons   
    names = get(findobj(gcf,'style','radiobutton'),'Tag');
    for i = 1:numel(names)
        pushed = sum( strcmp( SET(:,2), names{i} ) );
        if pushed
            handles.(names{i}).Value = 1;
            if isempty(findstr(names{i},'FromFile'))
                func( [handles.(names{i}).Tag '_Callback'], handles.(names{i}),'load',handles)
            end
        end
    end
    
    %% Set values in toggle buttons   
    names = get(findobj(gcf,'style','togglebutton'),'Tag');
    for i = 1:numel(names)
        pushed = sum( strcmp( SET(:,2), names{i} ) );
        if pushed
            handles.(names{i}).Value = 1;
            func( [handles.(names{i}).Tag '_Callback'], handles.(names{i}),'load',handles)
        end
    end
       
    %% Print
    fprintf(' Done!\n');
end