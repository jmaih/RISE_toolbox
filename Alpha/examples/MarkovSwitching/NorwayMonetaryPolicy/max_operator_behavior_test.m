function dev=max_operator_behavior_test(obj)
% this function checks that expected monetary policy shocks are positive
% whenever a the interest rate is expected to be below a cutoff value of
% 0.5

persistent targets shock_start shock_end mp_id

if isempty(targets)
    cutoff_interest_rate=-0.01;
    
    % 1- start and end date for the conditions
    start_date_cond=obj.options.cond_data_ct.TimeInfo(1);
    end_date_cond=obj.options.cond_data_ct.TimeInfo(end);
    % 2- start and end date for the smoothed
    start_date_shocks=obj.Filters.smoothed_shocks.regime_1.TimeInfo(1);
    end_date_shocks=obj.Filters.smoothed_shocks.regime_1.TimeInfo(end);
    % 3- overall start and end dates
    start_date=max(start_date_cond,start_date_shocks);
    end_date=min(end_date_cond,end_date_shocks);
    
    % location in each database
    cond_start=start_date_cond.date_2_observation(start_date);
    cond_end=start_date_cond.date_2_observation(end_date);
    
    % location in each database
    shock_start=start_date_shocks.date_2_observation(start_date);
    shock_end=start_date_shocks.date_2_observation(end_date);
    
    % search for the location of the interest rate
    R_id=find(strcmp('RN3M_NW',obj.options.cond_data_ct.varnames));
    % collect the corresponding conditioning information, which includes all
    % the pages of the database
    conditions=squeeze(cell2mat(obj.options.cond_data_ct.data(1+(cond_start:cond_end),R_id+1,:)));
    
    % locate the monetary policy shock in the list of the exogenous, which
    % should be the same as the order in the database
    mp_id=find(strcmp('EI',{obj.varexo.name}));
    
    % if we have used the iwb hypothesis, both conditions and monetary_policy_shocks
    % should have the same number of pages. Having matched the dates, we still
    % have to avoid possible nans
    targets=~isnan(conditions) & conditions<=cutoff_interest_rate;
% else
%     disp('persistent is effective')
end
% collect monetary policy shocks (current and expected).
monetary_policy_shocks=squeeze(cell2mat(obj.Filters.smoothed_shocks.regime_1.data(1+(shock_start:shock_end),mp_id+1,:)));

% now locate the expected monetary policy shocks that do not comply
dev=monetary_policy_shocks(targets)<0;

% % % % % rest_start=obj.options.restrictions.rest_start;
% % % % % % penalize the violations
% % % % % var_id=[obj.options.restrictions.restricted_variables.id];
% % % % % value=[obj.options.restrictions.restricted_variables.value];
% % % % % % enforce the sign here even if this is already taken care of
% % % % % % in set_options
% % % % % value=sign(value);
% % % % % % check how stephane does this
% % % % % additional_penalty=100;
% % % % % for ii=1:numel(var_id)
% % % % %     violations=sum(sum(sum(...
% % % % %         sign(Filters.eta(var_id(ii),:,rest_start:end))==value(ii)...
% % % % %         )));
% % % % %     logpost=logpost+violations*additional_penalty;
% % % % % end
