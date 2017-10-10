function self=recreate_parameters(self,lag_start)

nd=max(0,self.nx-self.constant);

p=abstvar.create_parameters_names(...
    self.ng*self.nvars,(lag_start:self.nlags),...
    self.ng*self.constant,self.ng*nd,self.markov_chains);

self.parameters=p.pnames(:).';

[self.mapping,self.markov_chains]=abstvar.map_estimation(p,...
    self.markov_chains,...
    self.endogenous,...
    self.exogenous,...
    self.members,self.nx);

do_sort=false;

[~,self.markov_chain_info]=parser.build_markov_regimes([],self.markov_chains,do_sort);

end