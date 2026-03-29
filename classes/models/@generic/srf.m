%--- help for generic/srf ---
%
%  SRF : shock response function, a flexible way to implement irfs.
% 
%  Syntax :
% 
%    myirfs=srf(m)
% 
%    myirfs=srf(m,x0)
% 
%    myirfs=srf(m,x0,r0)
% 
%    myirfs=srf(m,x0,r0,shkCombo)
% 
%    myirfs=srf(m,x0,r0,shkCombo,irf_periods)
% 
%    myirfs=srf(m,x0,r0,shkCombo,irf_periods,order)
% 
%    myirfs=srf(m,x0,r0,shkCombo,irf_periods,order,girf_runs)
% 
%  Inputs:
% 
%  - **m** [rise|dsge]: model object(s)
% 
%  - **x0** [empty|struct]: initial conditions for endogenous variables
% 
%  - **r0** [empty|integer]: initial/historical regime
% 
%  - **shkCombo** [numeric|struct|empty]: shock information such that 
% 
%    - if numeric : shock sign and or magnitude
% 
%    - if struct : structure with the names of the shocks and their values
%      for  the initial impulse : the impulse response is computed for a
%      **cocktail of shocks (scenario)** rather than for one specific shock.
%      Note that in this case, the outputwill be the same for all individual
%      shocks!!! 
% 
%  - **irf_periods** [empty|numeric]: period length for the irfs
% 
%  - **order** [empty|integer]: order of the solution used in the
%    computation of the irfs
% 
%  - **girf_runs** [empty|numeric]: if girf_runs>1, this triggers the
%    computation of generalized impulse responses. It automatically sets
%    **irf_regime_specific** and **irf_girf_regime_uncertainty** to false
% 
%  Outputs:
% 
%  - **myirfs** [struct]: structure in which the fields are the names of the
%    shocks of the model. This is true whether we have a cocktail of shocks
%    or not. The only difference is that for convenience, under a cocktail
%    of shocks, the responses are the same for all the shocks in the model.
%