function [ep]=format_transition_probabilities(ep,cn,nstates,vList,pList)

if isempty(ep)
    
    ep=constant_transitions(cn,nstates,pList);
    
    return
    
end

alltoks=union(vList,pList);

nv=numel(vList);

np=numel(pList);

v=rand(1,nv);  p=rand(1,np);

pat=['\<',parser.cell2matize(alltoks(:).'),'\>'];

probnames=abstvar.problist(cn,nstates);

if isempty(ep)
    
    return
    
end

ping=regexp(ep,'=','split');

nping=numel(ping);

nprobs=numel(probnames);

if nping~=nprobs
    
    error(sprintf(['expected %0.0f probabilities in ',...
        'markov chain %s but found %0.0f'],...
        nprobs,cn,numel(ping))) %#ok<SPERR>
    
end

myreplace=@do_replace;%#ok<NASGU> 
% @(x)num2str(rand);

for kk=1:nping
    
    left=ping{kk}{1};
    
    if ~ismember(left,probnames)
        
        error(sprintf(['%s is not a valid transition ',...
            'probability name for markov chain "%s"'],...
            left,cn)) %#ok<SPERR>
        
    end
    % remove the semicolon if any
    [right,ok]=encode_expression(ping{kk}{2}); % strrep(ping{kk}{2},';','');
    % add a semicolon to suppress output at evaluation
    if ~ok
        
        error(['"',ep{kk},'" is an invalid expression'])
        
    end
    
    ping{kk}{2}=right;
    
end

ep=[ping{:}];

tpn=ep(1:2:end);

tpv=cell2mat(strcat(ep(2:2:end),','));

tpv=str2func(['@(p,v)[',tpv(1:end-1),']']);

ep={tpn,tpv};

    function [neweqtn,ok]=encode_expression(eqtn)
        
        neweqtn=regexprep(eqtn,pat,'${myreplace($0)}');
        
        ok=is_valid_expression(neweqtn,p,v);
        
    end

    function newAtom=do_replace(atom)
        
        pORv='p';
        
        pos=find(strcmp(atom,pList));
        
        if isempty(pos)
            
            pORv='v';
            
            pos=find(strcmp(atom,vList));
            
        end
        
        newAtom=sprintf('%s(%0.0f)',pORv,pos);
        
    end

end

function  ep=constant_transitions(cn,nstates,pList)
    
tpn=abstvar.problist(cn,nstates);

tpv=num2cell(locate_variables(tpn,pList));

tpv=cellfun(@(x)sprintf('p(%0.0f)',x),tpv,'uniformOutput',false);

tpv=cell2mat(strcat(tpv(:).',','));

tpv=str2func(['@(p,v)[',tpv(1:end-1),']']);

ep={tpn,tpv};

end

function ok=is_valid_expression(express,p,v) %#ok<INUSD>

ok=true;

try
    
    result=evalc(express); %#ok<NASGU>
    
catch
    
    ok=false;
    
end

end
