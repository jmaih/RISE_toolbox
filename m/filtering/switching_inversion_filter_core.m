%  `switching_inversion_filter_core` implements a switching inversion filter for
%  regime-switching state space models. 
% 
%  The assume model is of the form
% 
%    y{t}=cy{r}+Z{r}*s{t}
% 
%    s{t}=cs{r}+T{r}*s{t-1}+R{r}*e{t}
% 
%    r is governed by transition Q(s{t})
% 
%    e{t} ~ N(0,I)
% 
%  INPUTS:
%  - `Q`: function handle computing the transition matrix.
%  - `T`: Transition matrix for the state space model (m,m,h).
%  - `R`: Covariance matrix for the state evolution noise (m,p,h).
%  - `Z`: Observation matrix mapping states to observations (p,m,h).
%  - `cs`: Constant term for state equation (m,h).
%  - `cy`: Constant term for observation equation (p,h).
%  - `y`: Observation data (p x n).
%  - `init`: Structure containint the Initial conditions for the filter.
% 
%  OUTPUTS:
%  - `loglik`: Log-likelihood of the data.
%  - `incr`: Log-likelihood increments at each time step.
%  - `sttm1r`: Filtered state vectors at t-1|t-1 for each regime.
%  - `sttr`: Filtered state vectors at t|t for each regime.
%  - `ettr`: Forecast errors for each regime.
% 
%  The `switching_inversion_filter_core` function performs filtering for state
%  space models with switching regimes.  It calculates the log-likelihood of
%  the data, log-likelihood increments at each time step, filtered state 
%  vectors at t-1|t-1, filtered state vectors at t|t, and forecast errors.
% 
%  The function iterates through each time step, updating state vectors and
%  probabilities for each regime. It uses the Kalman filter to estimate
%  states and log-likelihood at each time step. 
% 
%  See also: aggregate_collapse.
%