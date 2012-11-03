function flag=rise_isa(string,type)
if iscell(string)
    flag=false(size(string));
    for ii=1:numel(string)
        flag(ii)=engine(string{ii});
    end
else
    flag=engine(string);
end

    function flag=engine(string)
        switch type
            case 'atom'
                flag=isempty(regexp(string,'[/*\-+^()]','start'));
            case 'number'
                flag=isempty(regexp(string,'[^\d^\.]','start','once'));
            case 'function'
                flag=false;
                [start,finish]=regexp(string,'[a-zA-Z]+\w*\_*','start','end','once');
                if ~isempty(start) && isequal(start,1)
                    string=string(finish+1:end);
                    flag=isempty(string)||(strcmp(string(1),'(') && strcmp(string(end),')'));
                end
            otherwise
        end
    end
end

