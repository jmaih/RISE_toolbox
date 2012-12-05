function RawFile=read_file(FileName)
% this would be the equivalent of dynare's macroprocessor
SPACE_DELIMITERS=char([9:13,32]);
tank_open=false;
tank=cell(0,3);
RawFile=cell(0,1);
fid = fopen(FileName);
iter=0;
inclusions=false;
while 1
    rawline = fgetl(fid);
    if ~ischar(rawline), break, end
    rawline=remove_comments(rawline);
    if all(isspace(rawline))
        continue
    end

    iter=iter+1;
    % check whether there are files included by parsing rawline
    [tokk,rest_]=strtok(rawline,SPACE_DELIMITERS);
    if tank_open
        if strcmp(tokk,'@#endfor')
            rawline=cell(0,3);
            for ii=start:finish
                tank0=tank;
                tank0(:,1)=strrep(tank0(:,1),index,int2str(ii));
                rawline=[rawline;tank0];
            end
            tank_open=false;
            tank=cell(0,3);
        else
            tank=[tank;{rawline,FileName,iter}]; %#ok<*AGROW>
            continue
        end
    elseif strcmp(tokk,'@#include')
        quotes=strfind(rest_,'"');
        if isempty(rest_)
            error([mfilename,':: missing quotes after statement @#include in ',FileName,' at line',int2str(iter)])
        end
        if numel(quotes)~=2
            error([mfilename,':: expecting a pair of quotes in ',FileName,' at line ',int2str(iter)])
        end
        newfile=strtrim(rest_(quotes(1)+1:quotes(2)-1));
        if ~exist(newfile,'file')
            error([mfilename,':: file ',newfile,' not found::  ',FileName,' at line ',int2str(iter)])
        end
        rawline=read_file(newfile);
        inclusions=true;
    elseif strcmp(tokk,'@#for')
        inclusions=true;
        [index,rest_]=strtok(rest_,SPACE_DELIMITERS);
        index=['@{',index,'}'];
        in=strfind(rest_,'in');
        rest_=rest_(in+2:end);
        colon=strfind(rest_,':');
        start=eval(rest_(1:colon-1));
        finish=eval(rest_(colon+1:end));
        tank_open=true;
        continue
    elseif strcmp(tokk,'@#endfor')
        error([mfilename,':: preprocessor failed ',FileName,' at line ',int2str(iter)])
    else
        rawline={rawline,FileName,iter};
    end
    RawFile=[RawFile;rawline];
end
fclose(fid);

if inclusions
    expanded=['expanded_',FileName(1:end-4),FileName(end-3:end)];
    try %#ok<TRYNC> % this would fail under parallel computing
        if exist(expanded,'file')
            delete(expanded)
        end
        fid=fopen(expanded,'w+');
        for ii=1:size(RawFile,1)
            fprintf(fid,'%s \n',RawFile{ii,1});
        end
        fclose(fid);
    end
end
end

function rawline_=remove_comments(rawline_)
% locate comments
loc_=strfind(rawline_,'%');
if ~isempty(loc_)
    rawline_=rawline_(1:loc_(1)-1);
end
loc_=strfind(rawline_,'//');
if ~isempty(loc_)
    rawline_=rawline_(1:loc_(1)-1);
end
end

