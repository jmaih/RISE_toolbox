%  Matlab complains that there is no documentation. But then here is one!!!
%  To do:
%  - beggar-thy-neighbor: whenever you sample a new guy, you put him in the
%  striation it belongs i.e. in the striation with the same level of energy
%  - striation specific metropolis, i.e. each striation has its own
%  metropolis point, initialized at each iteration to be the mean of the
%  guys in the striation: this breaks correlation
%  - re-introduce geometric weights such that the guys with low energy have
%  a wider striation?
%  - we could also consider randomizing the weights or having weights depend
%  on the performance?
%  - finish the equi-energy sampler, doing it with structures like in this
%  function. check the adaptive equi-energy sampler for tips for
%  improvement.
%  - finish the parallel tempering sampler and optimize it, again write in
%  as structures. check for the latest developments
%  - re-write the metropolis hastings algorithm with the same philosophy as
%  in this function and make it independent of rise objects
%  - re-write the marginal data densities: 
%    - Meng and Wong
%    - Adaptive Importance Sampling,
%    - Mueller, 
%    - Chib-Jeliazkov, 
%    - Waggoner and Zha, 
%    - Geweke
%  - Include moves from the metaheuristic literature in the sampler? I most
%  likely would be blamed for insulting convergence theorems and
%  reversibility of markov processes?
%  - Do a paper that
%    - presents the problem
%    - scale intervals if the product goes to infinity?
%    - argues that the techniques are also and perhaps more useful for
%    nonlinear models estimated using frequentist approaches
%    - runs examples on two models
%        - the Schorfheide model
%        - the Liu-Waggoner-Zha model. Try to see whether the simulator
%        recover the posterior mode
%    - plot both the bivariate marginal and univariate marginal and argue
%    that they could be multimodal, while the metropolis would want to have
%    them unimodal.
%  - Remove field posterior simulation ? where to put the posterior
%  statistics?
%  - Does it have to be the case that the number of draws should hold in
%  memory for all the exercises including the calculation of MDD, IRFS,
%  intervals, etc. or does it have to be like dynare where we collect the
%  info from the disk as necessary?
%  - simulators independent
%  - put metropolis in the simulators department alongside EE, DSMH, etc.
%  - parfor with condition all intensive processes:
%    - first detect the presence of slaves
%    - curvature
%    - derivatives?
%    - plots
%    - scale 
%    - curvature
%