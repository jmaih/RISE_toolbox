function write_function_to_disk(fname,equations)

if isempty(equations)
    return
end

if isstruct(equations)
    neq=numel(equations);
    outstring='res1';
    strings='';
    closing_strings='end';
    for ii=1:neq
        if ii>1
            outstring=[outstring,',res',int2str(ii)]; %#ok<AGROW>
            closing_strings=char(closing_strings,'end');
        end
        strings=char(strings,['if nargout >',int2str(ii-1)]);
        strings=char(strings,equations(ii).formatted{1});
        strings=char(strings,['res',int2str(ii),'=res;']);
    end
    strings=char(strings,closing_strings);
else
    % if equations contains individual elements with more than one entry, we
    % need to transform and re-transform
    if size(equations{1},1)>1
        equations=char(equations);
        equations=cellstr(equations);
    end
    neq=numel(equations);
    newversion=false;
    if newversion
        strings='res=[';
        for eq=1:neq
            strings=char(strings,equations{eq});
        end
        strings=char(strings,'];');
    else
        disp([mfilename,...
            ':: this way of writing functions does not work with automatic derivatives'])
        strings='k=size(y,2);';
        strings=char(strings,['res= nan(',int2str(neq),',k);']);
        for eq=1:neq
            strings=char(strings,['res(',int2str(eq),',:)=',equations{eq},';']);
        end
    end
    outstring='res';
end

fid=fopen([fname,'.m'],'w');
fprintf(fid,'%s \n\n',['function [',outstring,']= ',fname,'(y,x,param,ss,def) %#ok<INUSD,INUSL>']);
for ii=1:size(strings,1)
    fprintf(fid,'%s \n\n',strings(ii,:));
end

fclose(fid);
