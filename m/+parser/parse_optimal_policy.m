function [dictionary,blocks,Model_block,PlannerObjective_block,jac_toc_]=...
    parse_optimal_policy(Model_block,dictionary,blocks)

current_block_id=find(strcmp('planner_objective',{blocks.name}));

% auxiliary variables and equations may also come from the planner
% objective. In that case we need to add auxiliary equations to the model
% block
[dictionary,PlannerObjective_block,Model_block]=...
    parser.planner_objective(dictionary,blocks(current_block_id).listing,...
    Model_block);

% remove item from block
blocks(current_block_id)=[];

jac_toc_=[];
dictionary.is_optimal_simple_rule_model=...
    ~dictionary.is_deficient_eqtns && dictionary.is_model_with_planner_objective;

dictionary.is_optimal_policy_model=...
    dictionary.is_model_with_planner_objective && ...
    dictionary.is_deficient_eqtns;

if dictionary.is_optimal_policy_model|| dictionary.is_optimal_simple_rule_model
    
    [Model_block,dictionary,jac_toc_]=parser.optimal_policy_system(...
        PlannerObjective_block,Model_block,dictionary);
    
end

end
