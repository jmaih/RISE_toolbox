function [dictionary,PlannerObjective_block,is_model_with_planner_objective]=planner_objective(dictionary,listing)
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


is_model_with_planner_objective=~isempty(listing);
PlannerObjective_block=cell(0,1);

dictionary.planner_system=[];
if ~is_model_with_planner_objective
    return
end

% there can't be multiple planner_objective blocks
% pull everything on one line
fichier_name_=listing{1,3};
line_number_=listing{1,1};
OneLiner='';
for ii=1:size(listing,1)
    OneLiner=[OneLiner,listing{ii,2}]; %#ok<*AGROW>
    if ~ismember(listing{ii,1},line_number_)
        line_number_=[line_number_,listing{ii,1}];
    end
end
OneLiner(isspace(OneLiner))=[];
Commitment='commitment=1;';
Discount='discount=.99;';
Objective='';
[startd,endd]=regexp(OneLiner,'discount=\w*\.?\w+','start','end');
if ~isempty(startd)
    Discount=[OneLiner(startd:endd),';'];
    OneLiner(startd:endd)=[];
end
[startc,endc]=regexp(OneLiner,'commitment=\w*\.?\w+','start','end');
if ~isempty(startc)
    Commitment=[OneLiner(startc:endc),';'];
    OneLiner(startc:endc)=[];
    % Time to decide whether to add a markov chain controlling policy
    % behavior or not. This depends on whether Commitment is loose
    equality=strfind(Commitment,'=');
    true_commitment=Commitment(equality+1:end-1);
    if ~ismember(true_commitment,{'0','1'})
        dictionary.markov_chains(end+1)=parser.initialize_markov_chain(...
            parser.loose_commit(),2,'is_endogenous',false);
        % Adding a new markov chain through loose commitment. In
        % re-ordering the chains, we need to update the info in the
        % parameters
        [~,tmp]=sort({dictionary.markov_chains.name});
        for iparam=1:numel(dictionary.parameters)
            dictionary.parameters(iparam).governing_chain=...
                find(dictionary.parameters(iparam).governing_chain==tmp);
        end
        dictionary.markov_chains=dictionary.markov_chains(tmp);
    end
end
if ~isempty(OneLiner)
    right_curly_brace=strfind(OneLiner,'}');
    if ~isempty(right_curly_brace)
        OneLiner=OneLiner(right_curly_brace+1:end);
    end
    Objective=OneLiner;
end
if strcmp(Objective(1),';')
    error([mfilename,':: No objective function provided. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
end
tmp={line_number_,Objective,fichier_name_; % policy objective
    line_number_,Commitment,fichier_name_; % commitment
    line_number_,Discount,fichier_name_};  % discount

[PlannerObjective_block,dictionary]=parser.capture_equations(dictionary,tmp,'planner_objective');

% now make sure there are variables and there is no time in the variables
time_vars=cell2mat(PlannerObjective_block{1,1}(2,~cellfun(@isempty,PlannerObjective_block{1,1}(2,:))));
if isempty(time_vars)
    error([mfilename,':: Planner objective must include variables. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
elseif any(time_vars)
    error([mfilename,':: leads or lags not allowed in planner_objective block. Check ',fichier_name_,' at line(s) ',mat2str(line_number_)])
end
dictionary.planner_system.model={Objective;Commitment;Discount};