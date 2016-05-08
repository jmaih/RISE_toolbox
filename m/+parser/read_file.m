function RawFile=read_file(FileName)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


RawFile=cell(0,1);

if isempty(FileName)
    
    return
    
end

fid = fopen(FileName);

iter=0;

is_block_comments_open=false;

while 1
    
    rawline = fgetl(fid);
    
    if ~ischar(rawline), break, end
    
    tokk=strtok(rawline);
    
    ibco=strcmp(tokk,'%{')||strcmp(tokk,'/*');
    
    if ibco 
        
        if is_block_comments_open
            
            error('nested block comments not allowed ')
            
        end
        
        is_block_comments_open=true;
        
    end
    
    if is_block_comments_open
        
        if strcmp(tokk,'%}')||strcmp(tokk,'*/')
            
            is_block_comments_open=false;
            
        end
        
        continue
        
    end
    
    rawline=parser.remove_comments(rawline);
    
    iter=iter+1;
    
    if all(isspace(rawline))
        
        continue
        
    end
    
    rawline={rawline,FileName,iter}; %#ok<*AGROW>
    
    RawFile=[RawFile;rawline];
    
end

fclose(fid);

end
