%  INTERNAL FUNCTION: creates impact matrix for contemporaneous and future shocks
% 
%  ::
% 
%    M=stochastic_impact_cumulator(model,y0,nsteps)
%    M=stochastic_impact_cumulator(model,y0,nsteps,y_pos)
%    M=stochastic_impact_cumulator(model,y0,nsteps,y_pos,e_pos)
%    M=stochastic_impact_cumulator(model,y0,nsteps,y_pos,e_pos,states)
%    [M,ufkst,states,PAI,TT,Q]=stochastic_impact_cumulator(...)
% 
%  Args:
% 
%     - **model** [struct]:
% 
%       - **T** [1 x h cell]: solution of the model,
%         T{regime}=[Ty,Tsig,Te_0,...Te_k]
%       - **sstate** [1 x h cell]: steady state in each regime
%       - **state_cols** [vector|{1:ny}]: location of state endogenous
%         variables in y0 (see below).
%       - **k** [scalar]: anticipation horizon (Beyond the current period)
%       - **Qfunc** [empty|function_handle]: endogenous transition matrix
% 
%     - **y0** [ny x 1 vector]: initial conditions
%     - **nsteps** [integer]: number of forecast steps
%     - **y_pos** [vector]: location of restricted endogenous variables.
%     - **e_pos** [vector]: location of restricted shocks
%     - **states** [empty|nsteps x 1 vector]: states visited for each forecast
%       step.
% 
%  Returns:
%     :
% 
%     - **M** [struct]:
% 
%       - **R** [matrix]: convoluted restrictions of shocks stemming from the
%         restrictions on endogenous variables
%       - **ufkst** [matrix]: unconditional forecasts for the restricted
%         endogenous variables (Excluding the initial conditions!!!).
%       - **const** [vector]: impact of the constant (steady state + risk) for
%         the restricted endogenous variables.
%       - **S** [matrix]: direct restrictions on shocks
%       - **nshocks** [integer]: number of shocks
%       - **ny** [integer]: number of endogenous variables
% 
%     - **ufkst** [ny x (nsteps+1) matrix]: Unconditional forecasts mean (with
%       initial condition at the beginning!!!)
%     - **states** [nsteps x 1 vector]: regimes visited for each forecast step.
%     - **PAI** [h x nsteps matrix]: choice probabilities of regimes for each
%       step
%     - **TT** [matrix]: convolution of autoregressive terms for the
%       restrictions
%     - **Q** [h x h x nsteps array]: Time series of transition matrices
% 
%