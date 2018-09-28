function RawFile=read_file(FileName,remove_comments)
% INTERNAL FUNCTION
%

if nargin<2
    
    remove_comments=[];
    
end

if isempty(remove_comments),remove_comments=true; end

RawFile=cell(0,1);

if isempty(FileName)
    
    return
    
end

FileName=char2struct(FileName);

% safeguard against the include files that sends out char instead of struct

for ifile=1:numel(FileName)
    
    RawFile_i=reading_engine(FileName(ifile),remove_comments);
    
    RawFile=[RawFile;RawFile_i];
    
end

    function st=char2struct(st)
        if ischar(st)
            tmp=struct();
            lastdot=find(st=='.',1,'last');
            tmp.ext=st(lastdot:end);
            tmp.fname=st(1:lastdot-1);
            st=tmp;
            % st=regexp(st,'(?<fname>\w+)(?<ext>\.\w+)','names');
            % this needs fixing for cases like ../Day1/baseline.rz that are
            % not properly parsed
        end
    end

end

function RawFile=reading_engine(FileName,remove_comments)

RawFile=cell(0,1);

fid = fopen([FileName.fname,FileName.ext]);

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
    
    iter=iter+1;
    
    if is_block_comments_open
        
        if strcmp(tokk,'%}')||strcmp(tokk,'*/')
            
            is_block_comments_open=false;
            
        end
        
        if remove_comments
            
            continue
            
        end
        
    end
    
    if remove_comments
        
        rawline=parser.remove_comments(rawline);
        
        if all(isspace(rawline))
            
            continue
            
        end
        
    end
    
    rawline={rawline,FileName.fname,iter}; %#ok<*AGROW>
    
    RawFile=[RawFile;rawline];
    
end

fclose(fid);

end
