function [nonlinear_restrictions,is_linear_restriction,derived_params]=...
    nonlinear_restrictions_engine(...
    param_names,regimes,chain_names,governing_chain,RestrictionsBlock)

isnonlin=@(x)any(x=='<')||any(x=='>');

isnonlin_not_allowed=@(x)~isempty(strfind(x,'=='));

n_restr=numel(RestrictionsBlock);

is_linear_restriction=true(1,n_restr);

[expr,replace,convert_the_guy]=regexp_setup2(param_names,...
    governing_chain,chain_names,regimes);

for irow=1:n_restr
    
    % change pname_chain_state into pname(chain,state)
    %-------------------------------------------------
    eqtn=parser.param_name_to_param_texname(RestrictionsBlock{irow},chain_names);
    
    if isnonlin_not_allowed(eqtn)
        
        disp(eqtn)
        
        error('== not allowed in nonlinear restrictions')
        
    elseif isnonlin(eqtn)
        % change pname or pname(chain,state) into M(i,j)
        RestrictionsBlock{irow}=regexprep(eqtn,expr,replace);
        
        is_linear_restriction(irow)=false;
        
    else
        
        RestrictionsBlock{irow}=process_definition(eqtn);
        
    end
    
end

derived_params=RestrictionsBlock(is_linear_restriction);

RestrictionsBlock(is_linear_restriction)=[];

nonlinear_restrictions=reprocess_nonlinear_restrictions(RestrictionsBlock(:).');

    function eqtn=process_definition(eqtn)
        
        equality=find(eqtn=='=');
        
        lhs=eqtn(1:equality-1);
        
        rhs=regexprep(eqtn(equality+1:end),expr,replace);
        
        cn=[];
        
        statepos=[];
        
        leftpar=strfind(lhs,'(');
        
        pname_lhs=lhs;
        
        if ~isempty(leftpar)
            
            rightpar=strfind(lhs,')');
            
            comma=strfind(lhs,',');
            
            cn=lhs(leftpar+1:comma-1);
            
            statepos=lhs(comma+1:rightpar-1);
            
            pname_lhs=lhs(1:leftpar-1);
            
        end
        
        [~,aloc,col]=convert_the_guy(pname_lhs,cn,statepos);
        
        eqtn={aloc,col(:).',str2func(['@(M)',rhs])};
        
    end

end

function [expr,replace,convert_the_guy]=regexp_setup2(param_names,...
    governing_chain,chain_names,regimes)
% parse expressions such as "pname", "pname(chain,state)"
% negative lookahead
pnames=[parser.cell2matize(param_names),'\>'];

c2='(\()?'; % group if exist AND capture

opt_cnames=[parser.cell2matize(chain_names-'const'),'?'];
c3=opt_cnames;

c4='(,)?'; % '(?:,)?'

opt_digits='(\d+)?';% capture if exist
c5=opt_digits;

c6='(\))?'; % group if exist but do not capture

expr=[pnames,c2,c3,c4,c5,c6];

replace='${convert_the_guy($1,$2,$3,$4,$5,$6)}';

convert_the_guy=@do_conversion;
    
    function [c,aloc,col]=do_conversion(pname,left_par,cn,comma,statepos,right_par)
        
        aloc=locate_variables(pname,param_names);
        
        if isempty(cn)
           
            cn='const';
            
            statepos='1';
        
        end
        
        if ~strcmp(chain_names(governing_chain(aloc)),cn)
            
            error(['parameter "',pname,'" is not controlled by markov chain "',cn,'"'])
        
        end
        
        c_id= strcmp(cn,chain_names);
        
        col=find(regimes(:,c_id)==str2double(statepos));
        
        if isempty(col)
            
            error(['wrong state number for parameter "',pname,'"'])
        
        end
        
        c=['M(',int2str(aloc),',',int2str(col(1)),')'];
        
        if isempty(left_par)
            % in case a right parenthesis was immediately followed by a
            % parameter name
            c=[c,comma,right_par];
            
        end
    
    end

end

function [nonlcon,nconst]=reprocess_nonlinear_restrictions(nonlcon)
% parameter restrictions
% for the moment,I only allow the matrix of parameters. the idea is that if
% the restrictions are violated, evaluation should fail and the steady
% state should not be computed. In general, one could think of allowing the
% steady state to enter this game, but I already have enough problems with
% defining the steady state in markov switching, etc.

% transform the nonlinear constraints. I would like to keep the
% flexibility of knowning what parameters enter the constraints and
% so I do not do this in format parameters
nconst=numel(nonlcon);

for iconstr=1:nconst
    
    % remove semicolon
    nonlcon{iconstr}=strrep(nonlcon{iconstr},';','');
    
    % now remove inequalities
    cutoff_type={'>=','<=','>','<','='};
    
    for itype=1:numel(cutoff_type)
        
        cutoff_locs=strfind(nonlcon{iconstr},cutoff_type{itype});
        
        if ~isempty(cutoff_locs)
            
            cutoff_type=cutoff_type{itype};
            
            break
            
        end
        
    end
    
    if ~isempty(cutoff_locs)
        
        span=length(cutoff_type);
        
        left=nonlcon{iconstr}(1:cutoff_locs-1);
        
        right=nonlcon{iconstr}(cutoff_locs+span:end);
        
        switch cutoff_type
            
            case '>='
                
                nonlcon{iconstr}=[right,'-(',left,')-eps;'];
                
            case '<='
                
                nonlcon{iconstr}=[left,'-(',right,')-eps;'];
                
            case '>'
                
                nonlcon{iconstr}=[right,'-(',left,');'];
                
            case '<'
                
                nonlcon{iconstr}=[left,'-(',right,');'];
                
            case '='
                
                nonlcon{iconstr}=['abs(',left,'-(',right,'))-eps;'];
                
        end
        
    end
    
end

nonlcon=cell2mat(nonlcon(:)');

nonlcon=nonlcon(1:end-1);

if ~isempty(nonlcon)
    
    nonlcon=str2func(['@(M)[',nonlcon,']']);
    
end

end