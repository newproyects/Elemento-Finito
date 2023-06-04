
[filename, pathname] = uigetfile('*.mat', 'Pick a mat-file');        
        figNumber = findobj('Tag','text1');       
        
        if(filename ~= 0)
            fullname = strcat(pathname, filename);
                       
            load(fullname);
            file_size = strfind(filename, '.');
            
            numbers = regexp(filename, '[0-9]');
            
            if (~(isempty(numbers)))
            file_size = min(numbers);
            end
            
            for i = 1:file_size - 1                
                names(i) = filename(i);           
            end           
            
            strcat(names,'');
            
            str= ...
                ['Data from disk has been loaded...'];
            
        else
            str= ...
                ['Data from disk has not been loaded...'];
        end
        set(figNumber,'String',str);
        










