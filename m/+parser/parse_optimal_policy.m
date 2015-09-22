function [dictionary,blocks,Model_block,PlannerObjective_block,jac_toc_]=...
    parse_optimal_policy(Model_block,dictionary,blocks)

current_block_id=find(strcmp('planner_objective',{blocks.name}));

[dictionary,PlannerObjective_block]=...
    parser.planner_objective(dictionary,blocks(current_block_id).listing);

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
