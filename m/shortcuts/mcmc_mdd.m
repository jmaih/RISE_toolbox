function varargout = mcmc_mdd(varargin)
% MCMC_MDD -- computes various types of log marginal data density
%
% Syntax
% -------
% ::
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb,ub)
%
%   [log_mdd,extras] = MCMC_MDD(theta_draws,lb,ub,options)
%
% Inputs
% -------
%
% - **theta_draws** [struct]: with fields "f" and "x". Each parameter is
% defined as a structure, which means that theta_draws is a vector of
% structures. "x" is the parameter vector and "f" is the NEGATIVE of the
% log posterior kernel evaluated at "x". In case "f" is instead the log
% posterior kernel itself, option **maximization** below has to be set to
% "true".
%
% - **lb** [empty|vector]: lower bound of the search space. Necessary only
% for the swz algorithm. Conveniently replaced with the lower bounds of
% theta_draws if empty.
%
% - **ub** [empty|vector]: upper bound of the search space. Necessary only
% for the swz algorithm. Conveniently replaced with the upper bounds of
% theta_draws if empty.
%
% - **options** [struct]: with possible fields
%   - **log_post_kern** [function handle]: function computing the log
%   posterior kernel for a given parameter vector
%   - **center_at_mean** [{false}|true]: if true, the distribution is
%   centered at the mean. Else, it is centered at the mode, which should be
%   the maximum of the log posterior kernel in theta_draws
%   - **algorithm** [{mhm}|swz|mueller|bridge|is|ris|cj]: 
%       - **mhm** is the modified harmonic mean
%       - **swz** is the Sims, Waggoner and Zha (2008) algorithm
%       - **mueller** is the unpublished Mueller algorithm (see Liu,
%       Waggoner and Zha 2011). 
%       - **bridge** is the Meng and Wong (1996) algorithm. 
%       - **is** is the Importance sampling algorithm. 
%       - **ris** is the reciprocal importance sampling algorithm. 
%       - **cj** is the Chib and Jeliazkov (2001) algorithm. 
%   - **L** [{[]}|integer]: number of IID draws
%   - **maximization** [{false}|true]: Informs the procedure about whether
%   we have a maximization or a minimization problem.
%   - **debug** [{false}|true]: print useful information during estimation.
%   - **mhm_tau** [{(.1:.1:.9)}|vector|scalar]: truncation probabilities
%   for the MHM algorithm 
%   - **swz_pvalue** [{90}|scalar]: scalar for the computation of the lower
%   bound in the SWZ algorithm
%   - **bridge_TolFun** [numeric|{sqrt(eps)}]: convergence criterion in the
%   BRIDGE algorithm
%
% Outputs
% --------
%
% - **log_mdd** [numeric]: log marginal data density
%
% - **extras** [empty|struct]: further output from specific algorithms
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

[varargout{1:nargout}]=utils.marginal_data_density.mcmc_mdd(varargin{:});