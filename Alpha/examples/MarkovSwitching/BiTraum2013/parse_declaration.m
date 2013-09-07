%classdef parse_declaration < handle
function list=parse_declaration(string,blocknames)
list=struct('name',{},'description',{},'start',{},'end',{},'class',{});
% classes are as follows
% 1 endogenous
% 2 exogenous
% 3 parameters
% 4 logvars
% 5 observables
tok='';
it=0;
iterlist=0;
descripOpen=false;
tokOpen=false;
while it<length(string)
    it=it+1;
    if isspace(string(it))||strcmp(string(it),',')||strcmp(string(it),';')
        assign()
    elseif strcmp(string(it),'"')
        if tokOpen
           error('starting a description right after a name')
        elseif descripOpen
            assign()
        else
            descripOpen=~descripOpen;
        end
    else
        if isempty(tok)
            start=it;
        end
        tok=[tok,string(it)];
    end
end

    function assign()
        if ~isempty(tok)
            if ismember(tok,blocknames)
                error([tok,' cannot be a declaration'])
            end
            if descripOpen
                list(iterlist).description=tok;
                descripOpen=false;
            else
                iterlist=iterlist+1;
                list(iterlist).name=tok;
                % by default, the description is the name
                list(iterlist).description=tok;
                tokOpen=false;
                finish=it-1;
                list(iterlist).start=start;
                list(iterlist).finish=finish;
            end
            tok='';
        end
    end
end