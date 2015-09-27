function obj=setup_nonlinear_restrictions(obj)
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

RestrictionsBlock=obj.options.estim_nonlinear_restrictions;
if isstruct(RestrictionsBlock)
    RestrictionsBlock=RestrictionsBlock.original;
end
if ~isempty(RestrictionsBlock)
    RestrictionsBlock=cellfun(@(x)x(~isspace(x)),...
        RestrictionsBlock,'uniformOutput',false);
end

n_restr=size(RestrictionsBlock,1);
isnonlin=@(x)any(x=='<')||any(x=='>');
isnonlin_not_allowed=@(x)~isempty(strfind(x,'=='));

param_names=obj.parameters.name;
governing_chain=obj.parameters.governing_chain;
chain_names=obj.markov_chains.small_markov_chain_info.chain_names;
regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
derived_parameters=[];

bend_it_like_a_chiquita();

RestrictionsBlock=transpose(RestrictionsBlock(:));

nonlinear_restrictions=reprocess_nonlinear_restrictions(RestrictionsBlock);
obj.number_of_restrictions=struct('auxiliary',sum(is_linear_restriction),...
    'linear',nan,...
    'nonlinear',sum(~is_linear_restriction));

obj=add_to_routines(obj,'derived_parameters',derived_parameters,...
    'nonlinear_restrictions',nonlinear_restrictions);

    function bend_it_like_a_chiquita()
        [expr,replace]=regexp_setup();
        convert_the_guy=@do_conversion; %#ok<NASGU>
        is_linear_restriction=true(1,n_restr);
        for irow=1:n_restr
            eqtn=RestrictionsBlock{irow};
            if isnonlin_not_allowed(eqtn)
                disp(eqtn)
                error('== not allowed in nonlinear restrictions')
            elseif isnonlin(eqtn)
                RestrictionsBlock{irow}=regexprep(eqtn,expr,replace);
                is_linear_restriction(irow)=false;
            else
                RestrictionsBlock{irow}=process_definition(eqtn);
            end
        end
        derived_parameters=RestrictionsBlock(is_linear_restriction);
        RestrictionsBlock(is_linear_restriction)=[];
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
            [~,aloc,col]=do_conversion(pname_lhs,cn,statepos);
            eqtn={aloc,col(:).',str2func(['@(M)',rhs])};
        end
        function [c,aloc,col]=do_conversion(pname,cn,statepos)
            aloc=locate_variables(pname,param_names);
            if isempty(cn)
                cn='const';
                statepos='1';
            end
            if ~strcmp(chain_names(governing_chain(aloc)),cn)
                error(['parameter "',pname,'" is not controled by markov chain "',cn,'"'])
            end
            c_id= strcmp(cn,chain_names);
            col=find(regimes(:,c_id)==str2double(statepos));
            if isempty(col)
                error(['wrong state number for parameter "',pname,'"'])
            end
            c=['M(',int2str(aloc),',',int2str(col(1)),')'];
        end
        function c2m=cell2matize(list)
            c2m=cell2mat(strcat(list,'|'));
            c2m=['(',c2m(1:end-1),')'];
        end
        function [expr,replace]=regexp_setup()
            capt_pnames=cell2matize(param_names);
            noleft_par='(?:\()?'; % group if exist but do not capture
            noright_par='(?:\))?'; % group if exist but do not capture
            opt_cnames=[cell2matize(chain_names-'const'),'?'];
            opt_comma=',?';
            opt_digits='(\d+)?';% capture if exist
            expr=[capt_pnames,noleft_par,opt_cnames,opt_comma,opt_digits,noright_par];
            replace='${convert_the_guy($1,$2,$3)}';
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
