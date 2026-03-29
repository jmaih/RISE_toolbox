%  mdd : Constructor for marginal data density objects
% 
%  Syntax::
% 
%        obj=mdd(theta_draws,log_post_kern,lb,ub)
% 
%        obj=mdd(theta_draws,log_post_kern,lb,ub,subset)
% 
%        obj=mdd(theta_draws,log_post_kern,lb,ub,subset,H)
% 
%        obj=mdd(theta_draws,log_post_kern,lb,ub,subset,H,maximization)
% 
%  Inputs
% 
%     - **theta_draws** [char|struct]: sampling drawss
%       => In case of a char, this is the location of the folder
%       containing the draws organized as described below
%       => In case of a structure, the fields are "f" and "x".
%       Each parameter vector is defined as a structure, which
%       means that theta_draws is a vector of structures. "x" is
%       the parameter vector and "f" is the NEGATIVE of the log
%       posterior kernel evaluated at "x". In case "f" is
%       instead the log posterior kernel itself, option
%       **maximization** below has to be set to "true".
% 
%     - **lb** [empty|vector]: lower bound of the search space.
%       Necessary only for the swz algorithm. Conveniently
%       replaced with the lower bounds of theta_draws if empty.
% 
%     - **ub** [empty|vector]: upper bound of the search space.
%       Necessary only for the swz algorithm. Conveniently
%       replaced with the upper bounds of theta_draws if empty.
% 
%     - **subset** (cell array|{empty}): When not empty, subset is a
%       1 x 2 cell array in which the first cell contains a
%       vector selecting the columns to retain in each chain and
%       the second column contains the chains retained. Any or
%       both of those cell array containts can be empty. Whenever
%       an entry is empty, all the information available is
%       selected. E.g. subsetting with dropping and trimming
%       mysubs={a:b:c,[1,3,5]}. In this example, the first
%       element selected is the one in position "a" and
%       thereafter every "b" element is selected until we reach
%       element in position "c". At the same time, we select
%       markov chains 1,3 and 5.
% 
%     - **H** [empty|matrix]: Hessian matrix for the parameters
%       Necessary only for the laplace algorithm. If left empty,
%       finite differences are used.
% 
%      - **maximization** [{false}|true|empty]: Informs the
%        procedure about whether we have a maximization or a
%        minimization problem.
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz
%
%    Documentation for mdd
%       doc mdd
%
%