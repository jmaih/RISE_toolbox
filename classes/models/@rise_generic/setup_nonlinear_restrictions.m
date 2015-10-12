function obj=setup_nonlinear_restrictions(obj)
% setup_nonlinear_restrictions - sets nonlinear restrictions
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
% - uses estim_nonlinear_restrictions, which should be a cell array. Each
% item of the array is a string of the form
%   - 'f(p1,p2,...,pn)>=h(p1,p2,...,pn)' 
%   - 'f(p1,p2,...,pn)>h(p1,p2,...,pn)' 
%   - 'f(p1,p2,...,pn)<=h(p1,p2,...,pn)' 
%   - 'f(p1,p2,...,pn)<h(p1,p2,...,pn)' 
%   - 'pj=h(p1,p2,...,pn)' 
%
% - In some cases, the explicit name for some parameter pj is not known in
% advance. In that case the name has to be formed explicitly as follows:
%   - pj=coef(eqtn,vbl,lag)
%   - pj=coef(eqtn,vbl,lag,chain,state)
%
% - In the statements above,
%   - eqtn [digits|variable name]
%   - vbl [digits|variable name]
%   - lag [digits]
%   - chain [char]
%   - state [digits]
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
endo_names=obj.endogenous.name;

bend_it_like_a_chiquita();

RestrictionsBlock=transpose(RestrictionsBlock(:));

nonlinear_restrictions=reprocess_nonlinear_restrictions(RestrictionsBlock);
obj.number_of_restrictions=struct('auxiliary',sum(is_linear_restriction),...
    'linear',nan,...
    'nonlinear',sum(~is_linear_restriction));

obj=add_to_routines(obj,'derived_parameters',derived_parameters,...
    'nonlinear_restrictions',nonlinear_restrictions);

    function bend_it_like_a_chiquita()
        [expr1,repl1,convcoef]=regexp_setup1(endo_names); %#ok<NASGU>
        [expr,replace,convert_the_guy]=regexp_setup2(param_names,...
            governing_chain,chain_names,regimes);
        is_linear_restriction=true(1,n_restr);
        for irow=1:n_restr
            % remove any coef(...), turning it into "pname" or
            % "pname(chain,state)"
            %---------------------------------------------------
            eqtn=regexprep(RestrictionsBlock{irow},expr1,repl1);
            % change pname_chain_state into pname(chain,state)
            %-------------------------------------------------
            eqtn=parser.param_name_to_param_texname(eqtn,...
                chain_names);
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
            [~,aloc,col]=convert_the_guy(pname_lhs,cn,statepos);
            eqtn={aloc,col(:).',str2func(['@(M)',rhs])};
        end
    end
end

function [expr,replace,convcoef]=regexp_setup1(endo_names)
% let nc = no capture
replace='${convcoef($1,$2,$3,$4,$5)}';
nc1='(?:coef\()';
eqtn_='(\w+)'; %1
nc2='(?:,)?';
vbl='(\w+)?'; %2
nc3='(?:,)?';
lag='(\d+)?'; %3
nc4='(?:,)?';
chain='(\w+)?'; %4
nc5='(?:,)?';
state='(\d+)?'; %5
nc6='(?:\))';
expr=[nc1,eqtn_,nc2,vbl,nc3,lag,nc4,chain,nc5,state,nc6];
convcoef=@coef_converter;
    function out=coef_converter(eqtn_,vbl,lag,chain,state)
        lag=str2double(lag);
        if isempty(vbl)
            out=eqtn_;
            if ~isvarname(out)
                error('nonlinear restriction badly specified')
            end
        else
            if isnan(lag)
                error('nonlinear restriction badly specified')
            end
            % don't go it all the way. Because then we create
            % pname_chain_state. What we actually want is
            % pname(chain,state).
                vv={eqtn_,vbl,lag};
            out=coef.create_parameter_name(vv,endo_names);
            if ~isempty(chain)
                if isempty(state)
                    error('nonlinear restriction badly specified')
                end
                out=[out,'(',chain,',',state,')'];
            end
        end
    end
end

function [expr,replace,convert_the_guy]=regexp_setup2(param_names,...
    governing_chain,chain_names,regimes)
% parse expressions such as "pname", "pname(chain,state)"
% negative lookahead
pnames=[cell2matize(param_names),'(?!\w+)'];
nc1='(?:\()?'; % group if exist but do not capture
opt_cnames=[cell2matize(chain_names-'const'),'?'];
nc2='(?:,)?';
nc3='(?:\))?'; % group if exist but do not capture
opt_digits='(\d+)?';% capture if exist
expr=[pnames,nc1,opt_cnames,nc2,opt_digits,nc3];
replace='${convert_the_guy($1,$2,$3)}';
convert_the_guy=@do_conversion;
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
