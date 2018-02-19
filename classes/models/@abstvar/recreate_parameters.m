function self=recreate_parameters(self,lag_start)

nd=max(0,self.nx-self.constant);

self.param_guide=abstvar.create_parameters_names(...
    self.ng*self.nvars,(lag_start:self.nlags),...
    self.ng*self.constant,self.ng*nd,self.nonvar_parameters);

pList=self.param_guide.pnames(:).';

vList=self.endogenous;

self.parameters=pList;

markov_chains=self.markov_chains;

self.is_time_varying_trans_prob=false;

for ic=1:numel(markov_chains)
    
    if ~self.is_time_varying_trans_prob
        
        self.is_time_varying_trans_prob=...
            ~isempty(markov_chains(ic).endogenous_probabilities);
    
    end
    
    [markov_chains(ic).endogenous_probabilities]=...
        abstvar.format_transition_probabilities(...
        markov_chains(ic).endogenous_probabilities,...
        markov_chains(ic).name,markov_chains(ic).number_of_states,...
        vList,pList);

end

[self.mapping,self.markov_chains]=abstvar.parameters_solution_mapper(...
    self.param_guide,...
    markov_chains,...
    self.endogenous,...
    self.exogenous,...
    self.members,self.nx);

do_sort=true;

[~,self.markov_chain_info]=parser.build_markov_regimes([],self.markov_chains,do_sort);

if do_sort
    % If everything is sorted we can all go home early
    %-------------------------------------------------
    pos=locate_variables(self.markov_chain_info.chain_names,{self.markov_chains.name});
    
    self.markov_chains=self.markov_chains(pos);
    
end

end