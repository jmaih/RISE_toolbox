function obj=setup_identification(obj,RestrictionsBlock)

param_names=obj.parameters.name;

par_nbr=sum(obj.parameters.number);
regimes=cell2mat(obj.markov_chains.small_markov_chain_info.regimes(2:end,2:end));
reg_nbr=size(regimes,1);
if isempty(obj.parameter_values)
    obj.parameter_values=nan(par_nbr,reg_nbr);
end
chain_names=obj.markov_chains.small_markov_chain_info.chain_names;

n_restr=size(RestrictionsBlock,1);
for ii=1:n_restr
    equation_i=RestrictionsBlock{ii};
    for rhs=1:size(equation_i,2)
        if ~isempty(equation_i{2,rhs})
            pname=equation_i{1,rhs};
            state=equation_i{2,rhs}{2};
            chain_=equation_i{2,rhs}{1};
            p_id=find(strcmp(pname,param_names));
            c_id= strcmp(chain_,chain_names);
            col=find(regimes(:,c_id)==state,1,'first');
            % Now we just need to replace the correct parameter location in
            % the matrix...
            equation_i{1,rhs}=['M(',sprintf('%0.0f',p_id),',',sprintf('%0.0f',col),')'];
        end
    end
    RestrictionsBlock{ii}=cell2mat(equation_i(1,:));
end
%---------------------------------------------
% parameter restrictions
% for the moment,I only allow the matrix of parameters. the idea is that if
% the restrictions are violated, evaluation should fail and the steady
% state should not be computed. In general, one could think of allowing the
% steady state to enter this game, but I already have enough problems with
% defining the steady state in markov switching, etc.
obj.parameter_restrictions=transpose(RestrictionsBlock(:));
tmp=['@(M)4*any([',cell2mat(transpose(RestrictionsBlock(:))),']==0)'];
% restrictions are violated if the function returns 4, which is the retcode
if ~isempty(obj.parameter_restrictions)
    tmp(isspace(tmp))=[];
    semcols=strfind(tmp,';');
    tmp(semcols(end))=[]; % remove last semicolon
% Get rid of definitions
    tmp=parser.substitute_definitions(tmp,obj.definitions.shadow_dynamic);
end

% add it to the handle of functions
obj.routines.parameter_restrictions=str2func(tmp);
%---------------------------------------------