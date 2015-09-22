function obj=setup_parameter_restrictions(obj,RestrictionsBlock)
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


param_names=obj.parameters.name;

par_nbr=sum(obj.parameters.number);
regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
reg_nbr=size(regimes,1);
if isempty(obj.parameter_values)
    obj.parameter_values=nan(par_nbr,reg_nbr);
end
chain_names=obj.markov_chains.small_markov_chain_info.chain_names;

% nonlinear restrictions will be checked
% linear restrictions will be assigned...
n_restr=size(RestrictionsBlock,1);
isnonlin=@(x)any(x=='<')||any(x=='>');
isnonlin_not_allowed=@(x)~isempty(strfind(x,'=='));
is_linear_restriction=true(1,n_restr);
linear_restrictions=cell(n_restr,3);
for ii=1:n_restr
    equation_i=RestrictionsBlock{ii};
    % check whether it is an equality or an inequality restriction
    %-------------------------------------------------------------
    full_eqtn=cell2mat(equation_i(1,:));
    full_eqtn(isspace(full_eqtn))=[];
    if isnonlin_not_allowed(full_eqtn)
        error('== not allowed in nonlinear restrictions')
    elseif isnonlin(full_eqtn)
        is_linear_restriction(ii)=false;
    end
    if is_linear_restriction(ii)
        equal_loc=find(full_eqtn=='=');
        rhs_name=full_eqtn(1:equal_loc-1);
    end
    nonlin_flag=~is_linear_restriction(ii);
    linear_restriction_lhs_in_the_box=false;
    for rhs=1:size(equation_i,2)
        if ~isempty(equation_i{2,rhs})
            pname=equation_i{1,rhs};
            state=equation_i{2,rhs}{2};
            chain_=equation_i{2,rhs}{1};
            p_id=find(strcmp(pname,param_names));
            c_id= strcmp(chain_,chain_names);
            % Now we just need to replace the correct parameter location in
            % the matrix...
            col=find(regimes(:,c_id)==state);
            if nonlin_flag||linear_restriction_lhs_in_the_box
                equation_i{1,rhs}=['M(',sprintf('%0.0f',p_id),...
                    ',',sprintf('%0.0f',col(1)),')'];
            else
                if ~strcmp(rhs_name,pname)
                    error(['unable to successfully parse the equality restriction ',full_eqtn])
                end
                if linear_restriction_lhs_in_the_box
                    error('in linear restrictions, there can only be one lhs variable')
                end
                linear_restriction_lhs_in_the_box=true;
                linear_restrictions(ii,:)={p_id,col(:).',[]};
            end
        end
    end
    full_eqtn=cell2mat(equation_i(1,:));
    if is_linear_restriction(ii)
        eq_loc=find(full_eqtn=='=');
        linear_restrictions{ii,3}=str2func(['@(M)',full_eqtn(eq_loc+1:end-1)]);
    else
        RestrictionsBlock{ii}=full_eqtn;
    end
end
derived_parameters=linear_restrictions(is_linear_restriction,:);

RestrictionsBlock=RestrictionsBlock(~is_linear_restriction);
%---------------------------------------------
RestrictionsBlock=transpose(RestrictionsBlock(:));

nonlinear_restrictions=reprocess_nonlinear_restrictions(RestrictionsBlock);
obj.number_of_restrictions=struct('auxiliary',sum(is_linear_restriction),...
    'linear',nan,...
    'nonlinear',sum(~is_linear_restriction));

obj=add_to_routines(obj,'derived_parameters',derived_parameters,...
    'nonlinear_restrictions',nonlinear_restrictions);

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
