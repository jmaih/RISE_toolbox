function strout=valid_varnames(names)

vnames={'x','p','ss','y','param','def'};
if nargin==0
    strout=vnames;
else
    if ischar(names)
        names=cellstr(names);
    elseif ~iscellstr(names)
        error([mfilename,':: input must be char or cellstr'])
    end
    for iname=1:numel(names)
        str=names{iname};
        flag=false;
        under_score=strfind(str,'_');
        if ~isempty(under_score) && numel(under_score)==1
            if any(strcmp(str(1:under_score-1),vnames)) % <--- ismember(str(1:under_score-1),{'x','y','p','ss'})
                d=str2double(str(under_score+1:end));
                if ~isnan(d) && isequal(d,floor(d))
                    flag=true;
                end
            end
        end
        if ~flag
            error([mfilename,':: ',str,' is not a valid variable name. Variable names can only be of the form (y|x|p|ss)_\d*'])
        end
    end
    strout=names;
end

