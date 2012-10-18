function flag=is_atom(string)
% atom:= does not contain (+-*/^,)
method=3;
forbidden='()+-*/^,';
switch method
    case 1
        flag=~any(ismember(string,forbidden)); % 0.502 s
    case 2
        string=unique(string);
        flag=true;
        for ii=1:length(string)
            aa=string(ii);
            for jj=1:length(forbidden)
                if strcmp(aa,forbidden(jj))
                    flag=false;
                    return
                end
            end
        end
    case 3
        flag=true;
        for ii=length(string)
            flag=flag && ~any(string(ii)==forbidden);
            if ~flag
                break
            end
        end
    otherwise
        error([mfilename,':: method not implemented'])
end
end
