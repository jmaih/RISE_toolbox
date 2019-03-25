%  PURPOSE: Temporal disaggregation using the Chow-Lin method
%  ------------------------------------------------------------
%  SYNTAX: res=chowlin(Ylow,Xhigh,s,aggreg_type,estim_method);
%  ------------------------------------------------------------
%  INPUT: Ylow: Nx1 ---> vector of low frequency data
%         Xhigh: nxp ---> matrix of high frequency indicators (without intercept)
%         aggreg_type: type of disaggregation
%             aggreg_type=1 or 'flow' ---> sum (flow)
%             aggreg_type=2 or 'average' ---> average (index)
%             aggreg_type=3 or 'last'  ---> last element (stock) ---> interpolation
%             aggreg_type=4 or 'first' ---> first element (stock) ---> interpolation
%         s: number of high frequency data points for each low frequency data points 
%         estim_method: estimation method: 
%             estim_method=0 ---> weighted least squares 
%             estim_method=1 ---> maximum likelihood with grid
%             estim_method=2 ---> maximum likelihood with fmincon
%  ------------------------------------------------------------
%  OUTPUT: res: a structure
%            res.meth    ='Chow-Lin';
%            res.aggreg_type      = type of disaggregation
%            res.estim_method    = method of estimation
%            res.N       = nobs. of low frequency data
%            res.n       = nobs. of high-frequency data
%            res.pred    = number of extrapolations
%            res.s       = frequency conversion between low and high freq.
%            res.p       = number of regressors (including intercept)
%            res.Ylow       = low frequency data
%            res.Xhigh       = high frequency indicators
%            res.y       = high frequency estimate
%            res.y_dt    = high frequency estimate: standard deviation
%            res.y_lo    = high frequency estimate: sd - sigma
%            res.y_up    = high frequency estimate: sd + sigma
%            res.u       = high frequency residuals
%            res.U       = low frequency residuals
%            res.beta    = estimated model parameters
%            res.beta_sd = estimated model parameters: standard deviation
%            res.beta_t  = estimated model parameters: t ratios
%            res.rho     = innovational parameter
%            res.aic     = Information criterion: AIC
%            res.bic     = Information criterion: BIC
%            res.val     = Objective function used by the estimation method
%            res.r       = grid of innovational parameters used by the estimation method
%  ------------------------------------------------------------
%  SEE ALSO: 
%  ------------------------------------------------------------
%