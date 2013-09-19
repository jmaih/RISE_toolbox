% RISE Toolbox
% Version Alpha_20130920 (R2012b) 20-Sep-2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toolbox        RISE
%
% Version        Alpha_20130920 20-Sep-2013
%
% Contents path  C:\Users\Junior\Documents\GitHub\RISE_toolbox\Alpha
%
%
%%%%%%%%%%%%%%%%%%%%   path:    %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help rise_startup">rise_startup</a> - else% gnu/linux
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+distributions   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+distributions\beta">classes\+distributions\beta</a>                     - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\cauchy">classes\+distributions\cauchy</a>                   - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\distribution_tests">classes\+distributions\distribution_tests</a>       - % Test of distributions
%   classes\+distributions\empirical_cdf            - (No help available)
%   classes\+distributions\find_bounds              - (No help available)
%   <a href="matlab:help classes\+distributions\find_hyperparameters">classes\+distributions\find_hyperparameters</a>     - happy_endings=[1 % FMINCON: First order optimality conditions satisfied.|| LSQNONLIN: converged to a solution.
%   <a href="matlab:help classes\+distributions\gamma">classes\+distributions\gamma</a>                    - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\hyperparameter_residuals">classes\+distributions\hyperparameter_residuals</a> - cdfn(pub,a,b,varargin{:})-(1-.5*alpha)]; % matlab_gamma_definition
%   <a href="matlab:help classes\+distributions\inv_gamma">classes\+distributions\inv_gamma</a>                - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\inv_wishart">classes\+distributions\inv_wishart</a>              - u=0; % <--- zeros(k,v)
%   <a href="matlab:help classes\+distributions\kernel_density">classes\+distributions\kernel_density</a>           -  smoothing density estimation
%   <a href="matlab:help classes\+distributions\laplace">classes\+distributions\laplace</a>                  - unction varargout=laplace(lowerquantileORmean,upperquantileORstdev,prob,c,d)% double exponential
%   <a href="matlab:help classes\+distributions\left_triang">classes\+distributions\left_triang</a>              - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\logistic">classes\+distributions\logistic</a>                 - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\lognormal">classes\+distributions\lognormal</a>                - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\normal">classes\+distributions\normal</a>                   - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\pareto">classes\+distributions\pareto</a>                   - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\right_triang">classes\+distributions\right_triang</a>             - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\truncated_normal">classes\+distributions\truncated_normal</a>         - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\uniform">classes\+distributions\uniform</a>                  - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\weibull">classes\+distributions\weibull</a>                  - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\wishart">classes\+distributions\wishart</a>                  - u=0; % <--- zeros(k,v)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+distributions\+unfinished   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+distributions\+unfinished\dirichlet">classes\+distributions\+unfinished\dirichlet</a>           - % one dirichlet at a time
%   <a href="matlab:help classes\+distributions\+unfinished\mv_normal">classes\+distributions\+unfinished\mv_normal</a>           - C=b; % assume this is the cholesky...
%   <a href="matlab:help classes\+distributions\+unfinished\truncated_mv_normal">classes\+distributions\+unfinished\truncated_mv_normal</a> - references:
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+mcf   %%%%%%%%%%%%%%%%%%%%
%
%   classes\+mcf\plot_correlation_patterns - (No help available)
%   <a href="matlab:help classes\+mcf\plot_correlation_scatter">classes\+mcf\plot_correlation_scatter</a>  - npar=numel(parameter_names);
%   <a href="matlab:help classes\+mcf\plot_smirnov">classes\+mcf\plot_smirnov</a>              - % density plots
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+msre_solvers   %%%%%%%%%%%%%%%%%%%%
%
%   classes\+msre_solvers\compute_preconditioner      - (No help available)
%   <a href="matlab:help classes\+msre_solvers\functional_iteration">classes\+msre_solvers\functional_iteration</a>        - based on the system X+inv(Aplus*X+A0)*Aminus=0. This
%   <a href="matlab:help classes\+msre_solvers\fwz_newton_system">classes\+msre_solvers\fwz_newton_system</a>           - We need to rewrite the system in the form
%   classes\+msre_solvers\integrate_structure         - (No help available)
%   <a href="matlab:help classes\+msre_solvers\newton_kronecker">classes\+msre_solvers\newton_kronecker</a>            - function [T,G0] = newton_kronecker(T0,Gplus01,A0,Aminus,Q,nn,h,frwz)
%   <a href="matlab:help classes\+msre_solvers\newton_kronecker_">classes\+msre_solvers\newton_kronecker_</a>           - alternative newton algorithm based on the expansion of the system
%   <a href="matlab:help classes\+msre_solvers\newton_kronecker_iteration">classes\+msre_solvers\newton_kronecker_iteration</a>  - This algorithm expands the system X+inv(Aplus*X+A0)*Aminus=0 as in
%   <a href="matlab:help classes\+msre_solvers\newton_kronecker_iteration_">classes\+msre_solvers\newton_kronecker_iteration_</a> - This algorithm expands the system X+inv(Aplus*X+A0)*Aminus=0 as in
%   <a href="matlab:help classes\+msre_solvers\newton_system">classes\+msre_solvers\newton_system</a>               - based on the derivatives of the system X+inv(Aplus*X+A0)*Aminus=0. This
%   <a href="matlab:help classes\+msre_solvers\newton_system_">classes\+msre_solvers\newton_system_</a>              - based on the derivatives of the system X+inv(Aplus*X+A0)*Aminus=0. This
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+msre_solvers\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+msre_solvers\private\A_times_X">classes\+msre_solvers\private\A_times_X</a> - computes A*X assuming the zero columns of A are deleted. the nonzero
%   <a href="matlab:help classes\+msre_solvers\private\X_times_A">classes\+msre_solvers\private\X_times_A</a> - computes X*A assuming the zero columns of A are deleted. the nonzero
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+obsolete   %%%%%%%%%%%%%%%%%%%%
%
%   classes\+obsolete\penalty_function - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+obsolete\@rise_anonymous   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+obsolete\@rise_anonymous\rise_anonymous">classes\+obsolete\@rise_anonymous\rise_anonymous</a> - % it is assumed that all objects in the vector share the same
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+obsolete\@sad_reverse   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+obsolete\@sad_reverse\char">classes\+obsolete\@sad_reverse\char</a>        - % char itself is already taken care of
%   classes\+obsolete\@sad_reverse\hessian     - (No help available)
%   <a href="matlab:help classes\+obsolete\@sad_reverse\jacobian">classes\+obsolete\@sad_reverse\jacobian</a>    - unction [Jac,... % jacobian in string form with auxiliary terms
%   <a href="matlab:help classes\+obsolete\@sad_reverse\sad_reverse">classes\+obsolete\@sad_reverse\sad_reverse</a> - % holds the list of signs such that if used on a specific object
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+obsolete\@sadiff   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+obsolete\@sadiff\char">classes\+obsolete\@sadiff\char</a>             - % search for 0 instead
%   <a href="matlab:help classes\+obsolete\@sadiff\diff">classes\+obsolete\@sadiff\diff</a>             - check that I do not need to output the object thanks to the
%   <a href="matlab:help classes\+obsolete\@sadiff\differentiate">classes\+obsolete\@sadiff\differentiate</a>    - this function
%   <a href="matlab:help classes\+obsolete\@sadiff\hessian">classes\+obsolete\@sadiff\hessian</a>          - build wrt right here so that things don't get confused with the ordering
%   <a href="matlab:help classes\+obsolete\@sadiff\jacobian">classes\+obsolete\@sadiff\jacobian</a>         - objectives is a function handle or an array of function handles
%   <a href="matlab:help classes\+obsolete\@sadiff\neat">classes\+obsolete\@sadiff\neat</a>             - str=[func,'(',varargin{:},')'];
%   <a href="matlab:help classes\+obsolete\@sadiff\print">classes\+obsolete\@sadiff\print</a>            - remove unused definitions
%   <a href="matlab:help classes\+obsolete\@sadiff\sadiff">classes\+obsolete\@sadiff\sadiff</a>           - lassdef sadiff % < handle
%   <a href="matlab:help classes\+obsolete\@sadiff\sadiff_test_____">classes\+obsolete\@sadiff\sadiff_test_____</a> - % choose level
%   <a href="matlab:help classes\+obsolete\@sadiff\setup">classes\+obsolete\@sadiff\setup</a>            - objectives is a function handle or an array of function handles
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+obsolete\@sadiff\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+obsolete\@sadiff\private\create_handle">classes\+obsolete\@sadiff\private\create_handle</a> - ndex=sprintf('%0.10g',operCount);
%   <a href="matlab:help classes\+obsolete\@sadiff\private\update_line">classes\+obsolete\@sadiff\private\update_line</a>   - f any(isletter(item))||strncmp(handle,indx,4) % do not waste time writing constant
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+ols   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+ols\ordinary_least_squares">classes\+ols\ordinary_least_squares</a> - ig2=RSS/T; % maximum likelihood estimator
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+parser   %%%%%%%%%%%%%%%%%%%%
%
%   classes\+parser\any2str                        - (No help available)
%   classes\+parser\append_file                    - (No help available)
%   <a href="matlab:help classes\+parser\build_markov_regimes">classes\+parser\build_markov_regimes</a>           - state_names=[state_names,[chain_names{ichain},'_',sprintf('%0.0f',irow)]]; %#ok<AGROW>
%   <a href="matlab:help classes\+parser\burry_probabilities">classes\+parser\burry_probabilities</a>            - pattern for forward-looking variables
%   <a href="matlab:help classes\+parser\capture_equations">classes\+parser\capture_equations</a>              - quation=initialize_equation();%cell(2,0);
%   <a href="matlab:help classes\+parser\capture_parameterization">classes\+parser\capture_parameterization</a>       - initialize output
%   <a href="matlab:help classes\+parser\check_markov_chains_adequacy">classes\+parser\check_markov_chains_adequacy</a>   -  that the markov chain is parameterized in all of its elements
%   classes\+parser\code2func                      - (No help available)
%   <a href="matlab:help classes\+parser\declarations2dictionary">classes\+parser\declarations2dictionary</a>        - UNTITLED11 Summary of this function goes here
%   <a href="matlab:help classes\+parser\delimiters">classes\+parser\delimiters</a>                     - UNTITLED4 Summary of this function goes here
%   <a href="matlab:help classes\+parser\determine_status">classes\+parser\determine_status</a>               - 'parameters','param' % same as definitions
%   <a href="matlab:help classes\+parser\endogenous">classes\+parser\endogenous</a>                     - format endogenous, parameters, observables, etc
%   <a href="matlab:help classes\+parser\exogenous">classes\+parser\exogenous</a>                      - format endogenous, parameters, observables, etc
%   <a href="matlab:help classes\+parser\file2blocks">classes\+parser\file2blocks</a>                    - UNTITLED3 Summary of this function goes here
%   <a href="matlab:help classes\+parser\greekify">classes\+parser\greekify</a>                       - so far I will just greekify names but later on, I might also greekify
%   <a href="matlab:help classes\+parser\initialize_dictionary">classes\+parser\initialize_dictionary</a>          - UNTITLED9 Summary of this function goes here
%   classes\+parser\is_transition_probability      - (No help available)
%   <a href="matlab:help classes\+parser\kron">classes\+parser\kron</a>                           - kronecker multiplication of cell arrays of strings
%   classes\+parser\listing                        - (No help available)
%   <a href="matlab:help classes\+parser\logvars2logvars">classes\+parser\logvars2logvars</a>                - old endo names will be useful if there is a problem in the model and has
%   <a href="matlab:help classes\+parser\look_around">classes\+parser\look_around</a>                    - UNTITLED5 Summary of this function goes here
%   <a href="matlab:help classes\+parser\observable">classes\+parser\observable</a>                     - format endogenous, parameters, observables, etc
%   <a href="matlab:help classes\+parser\param_name_to_valid_param_name">classes\+parser\param_name_to_valid_param_name</a> - change the parameter names from name(chain,state) to name_chain_state
%   <a href="matlab:help classes\+parser\parameter">classes\+parser\parameter</a>                      - format endogenous, parameters, observables, etc
%   <a href="matlab:help classes\+parser\parameters">classes\+parser\parameters</a>                     - format endogenous, parameters, observables, etc
%   <a href="matlab:help classes\+parser\parse">classes\+parser\parse</a>                          - by the way, I can still declare exogenous and make them observable at the
%   <a href="matlab:help classes\+parser\planner_objective">classes\+parser\planner_objective</a>              - there can't be multiple planner_objective blocks
%   <a href="matlab:help classes\+parser\preparse">classes\+parser\preparse</a>                       - int2str(x)=sprintf('%.0f',x)
%   classes\+parser\push_if_validated              - (No help available)
%   <a href="matlab:help classes\+parser\read_file">classes\+parser\read_file</a>                      - rawline={rawline,FileName,iter}; %#ok<*AGROW>
%   <a href="matlab:help classes\+parser\remove_comments">classes\+parser\remove_comments</a>                - locate comments
%   <a href="matlab:help classes\+parser\remove_definitions">classes\+parser\remove_definitions</a>             - est=regexp(tank(:,1),'@{[\w/*\-+^]+}','match');%@{\w+}
%   <a href="matlab:help classes\+parser\replace_definitions">classes\+parser\replace_definitions</a>            - % the parentheses here are not very efficient
%   classes\+parser\replace_steady_state_call      - (No help available)
%   <a href="matlab:help classes\+parser\setfield">classes\+parser\setfield</a>                       - UNTITLED10 Summary of this function goes here
%   <a href="matlab:help classes\+parser\string_mult">classes\+parser\string_mult</a>                    - check whether there is a semicolon at the end of any of the strings
%   <a href="matlab:help classes\+parser\transition_probabilities">classes\+parser\transition_probabilities</a>       - remove semicolon
%   <a href="matlab:help classes\+parser\valid_names_in_text">classes\+parser\valid_names_in_text</a>            - change the parameter names from name(chain,state) to name_chain_state
%   <a href="matlab:help classes\+parser\valid_param_name_to_tex_name">classes\+parser\valid_param_name_to_tex_name</a>   - change the parameter names from name_chain_state to name(chain,state) and
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+quasi_monte_carlo   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+quasi_monte_carlo\halton">classes\+quasi_monte_carlo\halton</a>          - Examples:
%   classes\+quasi_monte_carlo\latin_hypercube - (No help available)
%   <a href="matlab:help classes\+quasi_monte_carlo\sobol">classes\+quasi_monte_carlo\sobol</a>           - Examples:
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+rise_moments   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+rise_moments\recursive_moments">classes\+rise_moments\recursive_moments</a> - {
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@hdmr   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@hdmr\hdmr">classes\@hdmr\hdmr</a> - % objective is either : f and theta or a function that will help
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@hdmr\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@hdmr\private\orthonormal_polynomial">classes\@hdmr\private\orthonormal_polynomial</a> - arning('off')%#ok<WNOFF> %,'MATLAB:divideByZero'
%   <a href="matlab:help classes\@hdmr\private\theta_to_x">classes\@hdmr\private\theta_to_x</a>             - normalize
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@planar   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@planar\abs                 - (No help available)
%   classes\@planar\acos                - (No help available)
%   classes\@planar\acosh               - (No help available)
%   classes\@planar\and                 - (No help available)
%   classes\@planar\asin                - (No help available)
%   classes\@planar\asinh               - (No help available)
%   classes\@planar\atan                - (No help available)
%   classes\@planar\atanh               - (No help available)
%   <a href="matlab:help classes\@planar\char">classes\@planar\char</a>                - % unary functions
%   classes\@planar\commute             - (No help available)
%   <a href="matlab:help classes\@planar\compose_derivatives">classes\@planar\compose_derivatives</a> - use element by element operations in case we have vectors ?
%   classes\@planar\cos                 - (No help available)
%   classes\@planar\cosh                - (No help available)
%   classes\@planar\cot                 - (No help available)
%   <a href="matlab:help classes\@planar\diff">classes\@planar\diff</a>                - % numbers/vectors and variables which are not part of differentiation
%   classes\@planar\eq                  - (No help available)
%   classes\@planar\erf                 - (No help available)
%   classes\@planar\exp                 - (No help available)
%   classes\@planar\ge                  - (No help available)
%   classes\@planar\gt                  - (No help available)
%   <a href="matlab:help classes\@planar\if_elseif">classes\@planar\if_elseif</a>           - check whether all the 'second' elements are zero
%   classes\@planar\if_then_else        - (No help available)
%   classes\@planar\initialize          - (No help available)
%   classes\@planar\is_one              - (No help available)
%   classes\@planar\is_zero             - (No help available)
%   classes\@planar\isnumeric           - (No help available)
%   classes\@planar\kron                - (No help available)
%   classes\@planar\le                  - (No help available)
%   <a href="matlab:help classes\@planar\load_varlist">classes\@planar\load_varlist</a>        - retrieves the list of all the variables
%   classes\@planar\log                 - (No help available)
%   classes\@planar\log10               - (No help available)
%   classes\@planar\lt                  - (No help available)
%   classes\@planar\max                 - (No help available)
%   classes\@planar\min                 - (No help available)
%   classes\@planar\minus               - (No help available)
%   classes\@planar\mpower              - (No help available)
%   classes\@planar\mrdivide            - (No help available)
%   classes\@planar\mtimes              - (No help available)
%   <a href="matlab:help classes\@planar\multinary_operation">classes\@planar\multinary_operation</a> - obj=commit(obj);
%   classes\@planar\ne                  - (No help available)
%   classes\@planar\normalcdf           - (No help available)
%   classes\@planar\normalpdf           - (No help available)
%   classes\@planar\or                  - (No help available)
%   <a href="matlab:help classes\@planar\planar">classes\@planar\planar</a>              - varlist={'a','b','c','x','d','e'}; wrt={'b','x','e'};
%   <a href="matlab:help classes\@planar\plus">classes\@planar\plus</a>                - % Simplify x+(-y) in x-y
%   classes\@planar\set                 - (No help available)
%   classes\@planar\sign                - (No help available)
%   classes\@planar\sin                 - (No help available)
%   classes\@planar\sinh                - (No help available)
%   classes\@planar\sqrt                - (No help available)
%   classes\@planar\tan                 - (No help available)
%   classes\@planar\tanh                - (No help available)
%   <a href="matlab:help classes\@planar\uminus">classes\@planar\uminus</a>              - % Simplify -(-x) in x
%   classes\@planar\uplus               - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@planar\+tmp   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@planar\+tmp\arguments                          - (No help available)
%   classes\@planar\+tmp\close_differentiation_session      - (No help available)
%   <a href="matlab:help classes\@planar\+tmp\commit">classes\@planar\+tmp\commit</a>                             - ref=['ref_',rise_sym_main_map.tag,'_',sprintf('%.0f',Count),'_']; % adding an underscore latter so that we can use strrep, which should be faster, rather than regexp
%   <a href="matlab:help classes\@planar\+tmp\differentiate">classes\@planar\+tmp\differentiate</a>                      - % multiplying letters with numbers give numbers. We need to take the
%   <a href="matlab:help classes\@planar\+tmp\equation2rise_sym">classes\@planar\+tmp\equation2rise_sym</a>                  - % re-create the function
%   <a href="matlab:help classes\@planar\+tmp\extend_differentiation_list">classes\@planar\+tmp\extend_differentiation_list</a>        - wrt{ii}=['zZzZzZz_',sprintf('%0.0f',iter)];
%   classes\@planar\+tmp\initialize_differentiation_session - (No help available)
%   <a href="matlab:help classes\@planar\+tmp\print">classes\@planar\+tmp\print</a>                              - min_ncount) % minimum number of occurrences at which substitution occurs
%   <a href="matlab:help classes\@planar\+tmp\push">classes\@planar\+tmp\push</a>                               - case 'get_derivative' %
%   <a href="matlab:help classes\@planar\+tmp\recompose">classes\@planar\+tmp\recompose</a>                          - % the effective number of calls might be lower if the argument does
%   <a href="matlab:help classes\@planar\+tmp\rise_sym">classes\@planar\+tmp\rise_sym</a>                           - if strcmp(class(func),'rise_sym')%#ok<STISA> %isa(func,'rise_sym')
%   <a href="matlab:help classes\@planar\+tmp\swap_references">classes\@planar\+tmp\swap_references</a>                    - if ~isempty(old_ref) && ~strncmp(old_ref,'G',1)
%   <a href="matlab:help classes\@planar\+tmp\system">classes\@planar\+tmp\system</a>                             - % before creating the function, check which variables actually enter and use those as arguments
%   <a href="matlab:help classes\@planar\+tmp\trim">classes\@planar\+tmp\trim</a>                               - if nargin<3
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise\assign_estimates">classes\@rise\assign_estimates</a>                       - this is the routine called during estimation for assigning new estimates
%   classes\@rise\check_derivatives                      - (No help available)
%   <a href="matlab:help classes\@rise\check_optimum">classes\@rise\check_optimum</a>                          - unction retcode=check_optimum(obj,varlist)%,r0,c0,N
%   <a href="matlab:help classes\@rise\compute_steady_state">classes\@rise\compute_steady_state</a>                   - In this file, there are two ways the steady state can be computed:
%   <a href="matlab:help classes\@rise\counterfactual">classes\@rise\counterfactual</a>                         - shocks_db is a rise_time_series object with alternative history for the shocks
%   <a href="matlab:help classes\@rise\create_estimation_blocks">classes\@rise\create_estimation_blocks</a>               - this function separates the parameters to estimate into blocks controled
%   classes\@rise\draw_parameter                         - (No help available)
%   <a href="matlab:help classes\@rise\estimate">classes\@rise\estimate</a>                               - varargin: data,optimset,optimizer,hessian_type,estim_start_from_mode,estim_parallel
%   <a href="matlab:help classes\@rise\filter">classes\@rise\filter</a>                                 - obj=msre_linear_filter();%markov_switching_kalman_filter();
%   <a href="matlab:help classes\@rise\forecast">classes\@rise\forecast</a>                               - unction [cond_fkst_db,cond_fkst_mean_db,uncond_fkst_db]=forecast(obj,varargin)% historical_db
%   <a href="matlab:help classes\@rise\forecast_real_time">classes\@rise\forecast_real_time</a>                     - see also plot_real_time
%   <a href="matlab:help classes\@rise\get">classes\@rise\get</a>                                    - lseif ismember(lower(PropertyName),{'par_vals','parameters'})% <--strcmpi(PropertyName,'par_vals')
%   <a href="matlab:help classes\@rise\historical_decomposition">classes\@rise\historical_decomposition</a>               - PURPOSE: Computes historical decompositions of a DSGE model
%   <a href="matlab:help classes\@rise\initialize_solution_or_structure">classes\@rise\initialize_solution_or_structure</a>       - 'm_x',... % formerly T
%   <a href="matlab:help classes\@rise\irf">classes\@rise\irf</a>                                    - 2- If the model is estimated, one may want to draw from the distribution,
%   <a href="matlab:help classes\@rise\is_stable_system">classes\@rise\is_stable_system</a>                       - this function checks that the solved system is
%   <a href="matlab:help classes\@rise\isnan">classes\@rise\isnan</a>                                  - params=[params,pname]; %#ok<*AGROW>
%   <a href="matlab:help classes\@rise\lazy_report">classes\@rise\lazy_report</a>                            - % maximum number of rows and cols per figure
%   <a href="matlab:help classes\@rise\load_parameters">classes\@rise\load_parameters</a>                        - This loads the parameters. This allows the user to quickly
%   <a href="matlab:help classes\@rise\log_marginal_data_density">classes\@rise\log_marginal_data_density</a>              - simulation_folder=obj.folders_paths.simulations;%[obj.options.results_folder,filesep,'simulations'];
%   <a href="matlab:help classes\@rise\log_posterior_kernel">classes\@rise\log_posterior_kernel</a>                   - estim_hyperparams=obj.estim_hyperparams;
%   <a href="matlab:help classes\@rise\log_prior_density">classes\@rise\log_prior_density</a>                      - % alternatives are
%   <a href="matlab:help classes\@rise\monte_carlo_filtering">classes\@rise\monte_carlo_filtering</a>                  - pval_cutoff=0.1;
%   <a href="matlab:help classes\@rise\posterior_marginal_and_prior_densities">classes\@rise\posterior_marginal_and_prior_densities</a> - =================
%   <a href="matlab:help classes\@rise\posterior_simulator">classes\@rise\posterior_simulator</a>                    - 'mcmc_diagcov_adjust_coef',1e-5,... % COVt=COV_{t-1}+mcmc_diagcov_adjust_coef*eye(npar)
%   <a href="matlab:help classes\@rise\print_estimates">classes\@rise\print_estimates</a>                        - params=strcat(param_names,' : ',param_vals); This will not preserve
%   <a href="matlab:help classes\@rise\print_estimation_results">classes\@rise\print_estimation_results</a>               - recision='%8.4f';
%   <a href="matlab:help classes\@rise\print_solution">classes\@rise\print_solution</a>                         - % get the location of the variables: can be model specific
%   <a href="matlab:help classes\@rise\prior_plots">classes\@rise\prior_plots</a>                            - aveUnderName0=cell(1,nobj);%
%   <a href="matlab:help classes\@rise\report">classes\@rise\report</a>                                 - rep_type ='varendo','varexo','varobs','parameters','solution'
%   <a href="matlab:help classes\@rise\rise">classes\@rise\rise</a>                                   - % Description
%   <a href="matlab:help classes\@rise\set">classes\@rise\set</a>                                    - fields=fieldnames(value);
%   <a href="matlab:help classes\@rise\set_options">classes\@rise\set_options</a>                            - function obj=set_options(obj,varargin)
%   classes\@rise\set_properties                         - (No help available)
%   <a href="matlab:help classes\@rise\simulate">classes\@rise\simulate</a>                               - UNTITLED6 Summary of this function goes here
%   <a href="matlab:help classes\@rise\simulation_diagnostics">classes\@rise\simulation_diagnostics</a>                 - first determine the number of chains and the number  of matrices in each
%   <a href="matlab:help classes\@rise\solve">classes\@rise\solve</a>                                  - structrual_matrices not stored into the object in case we need to "get"
%   <a href="matlab:help classes\@rise\solve_alternatives">classes\@rise\solve_alternatives</a>                     - syntax is allobj=solve_alternatives(obj,solver,file2save2)
%   <a href="matlab:help classes\@rise\stoch_simul">classes\@rise\stoch_simul</a>                            - unction [oo_]=stoch_simul(obj,var_list,varargin)%,omega
%   <a href="matlab:help classes\@rise\svar_solve">classes\@rise\svar_solve</a>                             - obj=rise(svar_struct)
%   <a href="matlab:help classes\@rise\theoretical_autocorrelations">classes\@rise\theoretical_autocorrelations</a>           - unction [A,retcode]=theoretical_autocorrelations(obj,varargin)%,resolve_flag
%   <a href="matlab:help classes\@rise\theoretical_autocovariances">classes\@rise\theoretical_autocovariances</a>            - unction [A,retcode]=theoretical_autocovariances(obj,varargin)%,resolve_flag
%   <a href="matlab:help classes\@rise\variance_decomposition">classes\@rise\variance_decomposition</a>                 - PURPOSE: Computes variance decompositions of a MSRE model
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise\private\compute_definitions">classes\@rise\private\compute_definitions</a>                   - evaluate definitions
%   <a href="matlab:help classes\@rise\private\concatenate">classes\@rise\private\concatenate</a>                           - precision='%8.6f';
%   <a href="matlab:help classes\@rise\private\draw_parameters">classes\@rise\private\draw_parameters</a>                       - this function draws a parameter and assigns it.
%   <a href="matlab:help classes\@rise\private\dsge_var_irf">classes\@rise\private\dsge_var_irf</a>                          - compute PHIb, SIGb and ZZi
%   <a href="matlab:help classes\@rise\private\first_order_solver">classes\@rise\private\first_order_solver</a>                    -  options are 0 or 'msre_klein' (for constant-parameter models)
%   <a href="matlab:help classes\@rise\private\format_parameters">classes\@rise\private\format_parameters</a>                     - baseline calibration
%   <a href="matlab:help classes\@rise\private\generate_starting_point">classes\@rise\private\generate_starting_point</a>               - this function attempts to reduce the number of rejection of randomly
%   <a href="matlab:help classes\@rise\private\latex_model_file">classes\@rise\private\latex_model_file</a>                      - ar_list = prepare_list(get(model,'par_list')); % bold
%   <a href="matlab:help classes\@rise\private\load_data">classes\@rise\private\load_data</a>                             - unction [obj,issue,retcode]=load_data(obj,varargin)%,estimation_flag
%   <a href="matlab:help classes\@rise\private\load_functions">classes\@rise\private\load_functions</a>                        - Get rid of definitions
%   <a href="matlab:help classes\@rise\private\load_mode">classes\@rise\private\load_mode</a>                             - reload the start values for the estimated parameters from a structure.
%   <a href="matlab:help classes\@rise\private\parameters_posterior_moments">classes\@rise\private\parameters_posterior_moments</a>          - approximated median calculated as the median of the medians
%   classes\@rise\private\potential_scale_reduction             - (No help available)
%   <a href="matlab:help classes\@rise\private\save_filters">classes\@rise\private\save_filters</a>                          - % re-order in terms of dates x nsteps x nvars
%   <a href="matlab:help classes\@rise\private\second_order_solver">classes\@rise\private\second_order_solver</a>                   - Ix=speye(endo_nbr);
%   <a href="matlab:help classes\@rise\private\setup_calibration">classes\@rise\private\setup_calibration</a>                     - (No help available)
%   <a href="matlab:help classes\@rise\private\setup_identification">classes\@rise\private\setup_identification</a>                  - % Now we just need to replace the correct parameter location in
%   <a href="matlab:help classes\@rise\private\setup_measurement_errors">classes\@rise\private\setup_measurement_errors</a>              - % pick only the endogenous observables
%   <a href="matlab:help classes\@rise\private\setup_priors">classes\@rise\private\setup_priors</a>                          - unction obj=setup_priors(obj,MyPriors,error_control)% is_switching=[ParameterInfo.is_switching];
%   <a href="matlab:help classes\@rise\private\simulation_engine">classes\@rise\private\simulation_engine</a>                     - simulate zeroth order
%   classes\@rise\private\store_probabilities                   - (No help available)
%   classes\@rise\private\substitute_definitions                - (No help available)
%   <a href="matlab:help classes\@rise\private\substitute_definitions_in_definitions">classes\@rise\private\substitute_definitions_in_definitions</a> - Get rid of definitions
%   <a href="matlab:help classes\@rise\private\svar_create">classes\@rise\private\svar_create</a>                           - This function creates a rise object of the svar type from a rational
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_endo_priors   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_endo_priors\rise_endo_priors">classes\@rise_endo_priors\rise_endo_priors</a> - Shat % zero-frequency spectral density
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_nad   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_nad\rise_nad">classes\@rise_nad\rise_nad</a> - % numerical automatic differentiation
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_pair   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_pair\rise_pair">classes\@rise_pair\rise_pair</a> - iter=0 % actual size
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_report   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_report\rise_report">classes\@rise_report\rise_report</a> - % add some space between records for lisibility and for
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_sym   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@rise_sym\abs                                - (No help available)
%   classes\@rise_sym\acos                               - (No help available)
%   classes\@rise_sym\acosh                              - (No help available)
%   classes\@rise_sym\and                                - (No help available)
%   classes\@rise_sym\arguments                          - (No help available)
%   classes\@rise_sym\asin                               - (No help available)
%   classes\@rise_sym\asinh                              - (No help available)
%   classes\@rise_sym\atan                               - (No help available)
%   classes\@rise_sym\atanh                              - (No help available)
%   classes\@rise_sym\close_differentiation_session      - (No help available)
%   <a href="matlab:help classes\@rise_sym\commit">classes\@rise_sym\commit</a>                             - ref=['ref_',rise_sym_main_map.tag,'_',sprintf('%.0f',Count),'_']; % adding an underscore latter so that we can use strrep, which should be faster, rather than regexp
%   classes\@rise_sym\commute                            - (No help available)
%   <a href="matlab:help classes\@rise_sym\compose_derivatives">classes\@rise_sym\compose_derivatives</a>                - t14 = t13/x.args{1}; % <-- creates problems when x.args{1}=0
%   classes\@rise_sym\cos                                - (No help available)
%   classes\@rise_sym\cosh                               - (No help available)
%   <a href="matlab:help classes\@rise_sym\diff">classes\@rise_sym\diff</a>                               - d=rise_sym_main_map.zero; % rise_sym.push('get_zero');%d=rise_sym(0);
%   <a href="matlab:help classes\@rise_sym\differentiate">classes\@rise_sym\differentiate</a>                      - % multiplying letters with numbers give numbers. We need to take the
%   classes\@rise_sym\eq                                 - (No help available)
%   <a href="matlab:help classes\@rise_sym\equation2rise_sym">classes\@rise_sym\equation2rise_sym</a>                  - % re-create the function
%   classes\@rise_sym\erf                                - (No help available)
%   classes\@rise_sym\exp                                - (No help available)
%   <a href="matlab:help classes\@rise_sym\extend_differentiation_list">classes\@rise_sym\extend_differentiation_list</a>        - wrt{ii}=['zZzZzZz_',sprintf('%0.0f',iter)];
%   classes\@rise_sym\ge                                 - (No help available)
%   classes\@rise_sym\get_incidence                      - (No help available)
%   classes\@rise_sym\gt                                 - (No help available)
%   <a href="matlab:help classes\@rise_sym\if_elseif">classes\@rise_sym\if_elseif</a>                          - check whether all the 'second' elements are zero
%   classes\@rise_sym\if_then_else                       - (No help available)
%   classes\@rise_sym\initialize_differentiation_session - (No help available)
%   classes\@rise_sym\is_one                             - (No help available)
%   classes\@rise_sym\is_zero                            - (No help available)
%   classes\@rise_sym\isnumeric                          - (No help available)
%   classes\@rise_sym\kron                               - (No help available)
%   classes\@rise_sym\le                                 - (No help available)
%   <a href="matlab:help classes\@rise_sym\load_varlist">classes\@rise_sym\load_varlist</a>                       - retrieves the list of all the variables
%   classes\@rise_sym\log                                - (No help available)
%   classes\@rise_sym\log10                              - (No help available)
%   classes\@rise_sym\lt                                 - (No help available)
%   classes\@rise_sym\max                                - (No help available)
%   classes\@rise_sym\min                                - (No help available)
%   classes\@rise_sym\minus                              - (No help available)
%   classes\@rise_sym\mpower                             - (No help available)
%   classes\@rise_sym\mrdivide                           - (No help available)
%   classes\@rise_sym\mtimes                             - (No help available)
%   <a href="matlab:help classes\@rise_sym\multinary_operation">classes\@rise_sym\multinary_operation</a>                - {
%   classes\@rise_sym\ne                                 - (No help available)
%   classes\@rise_sym\normalcdf                          - (No help available)
%   classes\@rise_sym\normalpdf                          - (No help available)
%   classes\@rise_sym\or                                 - (No help available)
%   <a href="matlab:help classes\@rise_sym\plus">classes\@rise_sym\plus</a>                               - % Simplify x+(-y) in x-y
%   <a href="matlab:help classes\@rise_sym\print">classes\@rise_sym\print</a>                              - min_ncount) % minimum number of occurrences at which substitution occurs
%   <a href="matlab:help classes\@rise_sym\push">classes\@rise_sym\push</a>                               - case 'get_derivative' %
%   <a href="matlab:help classes\@rise_sym\recompose">classes\@rise_sym\recompose</a>                          - % the effective number of calls might be lower if the argument does
%   <a href="matlab:help classes\@rise_sym\rise_sym">classes\@rise_sym\rise_sym</a>                           - if strcmp(class(func),'rise_sym')%#ok<STISA> %isa(func,'rise_sym')
%   classes\@rise_sym\sign                               - (No help available)
%   classes\@rise_sym\sin                                - (No help available)
%   classes\@rise_sym\sinh                               - (No help available)
%   classes\@rise_sym\sqrt                               - (No help available)
%   <a href="matlab:help classes\@rise_sym\swap_references">classes\@rise_sym\swap_references</a>                    - if ~isempty(old_ref) && ~strncmp(old_ref,'G',1)
%   <a href="matlab:help classes\@rise_sym\system">classes\@rise_sym\system</a>                             - % before creating the function, check which variables actually enter and use those as arguments
%   classes\@rise_sym\tan                                - (No help available)
%   classes\@rise_sym\tanh                               - (No help available)
%   <a href="matlab:help classes\@rise_sym\trim">classes\@rise_sym\trim</a>                               - if nargin<3
%   <a href="matlab:help classes\@rise_sym\uminus">classes\@rise_sym\uminus</a>                             - % Simplify -(-x) in x
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_sym\+junk   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_sym\+junk\compact_form">classes\@rise_sym\+junk\compact_form</a> - % handle will push this information automatically.
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_time_series   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_time_series\addpages">classes\@rise_time_series\addpages</a>                  - % first batch
%   classes\@rise_time_series\aggregate                 - (No help available)
%   classes\@rise_time_series\allmean                   - (No help available)
%   classes\@rise_time_series\and                       - (No help available)
%   <a href="matlab:help classes\@rise_time_series\ar">classes\@rise_time_series\ar</a>                        - - B: vector of regression coefficients in the linear model Y = X*B.
%   <a href="matlab:help classes\@rise_time_series\automatic_model_selection">classes\@rise_time_series\automatic_model_selection</a> - unction FinalResults=automatic_model_selection(tsdata,endo_name,exo_names,options)%[final_mod,x]
%   <a href="matlab:help classes\@rise_time_series\bar">classes\@rise_time_series\bar</a>                       - et(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %...
%   classes\@rise_time_series\bsxfun                    - (No help available)
%   <a href="matlab:help classes\@rise_time_series\collect">classes\@rise_time_series\collect</a>                   - % if there are many names then check that nargin==1
%   classes\@rise_time_series\corrcoef                  - (No help available)
%   classes\@rise_time_series\cov                       - (No help available)
%   <a href="matlab:help classes\@rise_time_series\double">classes\@rise_time_series\double</a>                    - cell2mat does not render the correct size when input is
%   classes\@rise_time_series\drop                      - (No help available)
%   classes\@rise_time_series\dummy                     - (No help available)
%   <a href="matlab:help classes\@rise_time_series\exp">classes\@rise_time_series\exp</a>                       - Here it does not make sense to have names any more. But
%   <a href="matlab:help classes\@rise_time_series\hist">classes\@rise_time_series\hist</a>                      - warning('second argument of hist is not numeric and was ignored')
%   classes\@rise_time_series\horzcat                   - (No help available)
%   <a href="matlab:help classes\@rise_time_series\hpfilter">classes\@rise_time_series\hpfilter</a>                  -  filters a collection of time series.
%   classes\@rise_time_series\interpolate               - (No help available)
%   classes\@rise_time_series\intersect                 - (No help available)
%   <a href="matlab:help classes\@rise_time_series\jbtest">classes\@rise_time_series\jbtest</a>                    - -H: H=0 means null hypothesis ("the data are normally
%   classes\@rise_time_series\kurtosis                  - (No help available)
%   <a href="matlab:help classes\@rise_time_series\line">classes\@rise_time_series\line</a>                      - 10             'yyyy'                   2000
%   classes\@rise_time_series\log                       - (No help available)
%   classes\@rise_time_series\max                       - (No help available)
%   classes\@rise_time_series\mean                      - (No help available)
%   classes\@rise_time_series\min                       - (No help available)
%   classes\@rise_time_series\minus                     - (No help available)
%   classes\@rise_time_series\mode                      - (No help available)
%   classes\@rise_time_series\mpower                    - (No help available)
%   classes\@rise_time_series\mrdivide                  - (No help available)
%   classes\@rise_time_series\mtimes                    - (No help available)
%   classes\@rise_time_series\nan                       - (No help available)
%   classes\@rise_time_series\ones                      - (No help available)
%   classes\@rise_time_series\pages2struct              - (No help available)
%   <a href="matlab:help classes\@rise_time_series\plot">classes\@rise_time_series\plot</a>                      - et(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %...
%   <a href="matlab:help classes\@rise_time_series\plot_separate">classes\@rise_time_series\plot_separate</a>             - obj=this.subsref(S); % <--- obj=this(vnames{id}); does not work
%   <a href="matlab:help classes\@rise_time_series\plot_window">classes\@rise_time_series\plot_window</a>               - % db(vnames{1}) would not work since we are inside a method of the
%   <a href="matlab:help classes\@rise_time_series\plotyy">classes\@rise_time_series\plotyy</a>                    - linkaxes(ax,'x')
%   <a href="matlab:help classes\@rise_time_series\plus">classes\@rise_time_series\plus</a>                      - % Here it does not make sense to have names any more. But
%   classes\@rise_time_series\rand                      - (No help available)
%   classes\@rise_time_series\randn                     - (No help available)
%   classes\@rise_time_series\range                     - (No help available)
%   <a href="matlab:help classes\@rise_time_series\rdivide">classes\@rise_time_series\rdivide</a>                   - just to make division robust
%   <a href="matlab:help classes\@rise_time_series\regress">classes\@rise_time_series\regress</a>                   - - B: vector of regression coefficients in the linear model Y = X*B.
%   classes\@rise_time_series\reset_start_date          - (No help available)
%   <a href="matlab:help classes\@rise_time_series\rise_time_series">classes\@rise_time_series\rise_time_series</a>          - % throw away trailing nan observations
%   classes\@rise_time_series\skewness                  - (No help available)
%   <a href="matlab:help classes\@rise_time_series\spectrum">classes\@rise_time_series\spectrum</a>                  - autocovariances
%   classes\@rise_time_series\std                       - (No help available)
%   <a href="matlab:help classes\@rise_time_series\step_dummy">classes\@rise_time_series\step_dummy</a>                - implements a step dummy in the time series
%   <a href="matlab:help classes\@rise_time_series\subsasgn">classes\@rise_time_series\subsasgn</a>                  - unction this=subsasgn(this,s,b)% subsasgn
%   <a href="matlab:help classes\@rise_time_series\subsref">classes\@rise_time_series\subsref</a>                   - %                     elseif  strcmp(s(1).type,'()')
%   classes\@rise_time_series\sum                       - (No help available)
%   classes\@rise_time_series\times                     - (No help available)
%   <a href="matlab:help classes\@rise_time_series\transform">classes\@rise_time_series\transform</a>                 - type is one of the following:
%   classes\@rise_time_series\uminus                    - (No help available)
%   classes\@rise_time_series\var                       - (No help available)
%   <a href="matlab:help classes\@rise_time_series\window">classes\@rise_time_series\window</a>                    - this=window(this,StartDate,EndDate)
%   classes\@rise_time_series\zeros                     - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_time_series\private   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@rise_time_series\private\CombineDates - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@sad_forward   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@sad_forward\hessian     - (No help available)
%   <a href="matlab:help classes\@sad_forward\jacobian">classes\@sad_forward\jacobian</a>    - objectives is a function handle or an array of function handles
%   <a href="matlab:help classes\@sad_forward\sad_forward">classes\@sad_forward\sad_forward</a> - lassdef sad_forward %< handle
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@sad_reverse   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@sad_reverse\hessian">classes\@sad_reverse\hessian</a>     - call the jacobian in non-vectorized form
%   <a href="matlab:help classes\@sad_reverse\jacobian">classes\@sad_reverse\jacobian</a>    - %
%   <a href="matlab:help classes\@sad_reverse\sad_reverse">classes\@sad_reverse\sad_reverse</a> - %------------------
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\log_marginal_data_density   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\log_marginal_data_density\chib_jeliazkov">classes\log_marginal_data_density\chib_jeliazkov</a>         - % now sample J parameter vectors from the proposal density. I don't think I
%   classes\log_marginal_data_density\modified_harmonic_mean - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\AcceleratedVsNonAcceleratedSolutions   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\AcceleratedVsNonAcceleratedSolutions\howto">examples\AcceleratedVsNonAcceleratedSolutions\howto</a>    - % housekeeping
%   <a href="matlab:help examples\AcceleratedVsNonAcceleratedSolutions\howto_sw">examples\AcceleratedVsNonAcceleratedSolutions\howto_sw</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\BoundsToMoments   %%%%%%%%%%%%%%%%%%%%
%
%   examples\BoundsToMoments\howto - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\ChangingParametersOutsideModelFile   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\ChangingParametersOutsideModelFile\howto">examples\ChangingParametersOutsideModelFile\howto</a> - ,...
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\DsgeVar   %%%%%%%%%%%%%%%%%%%%
%
%   examples\DsgeVar\datarabanal_hybrid - (No help available)
%   <a href="matlab:help examples\DsgeVar\howto">examples\DsgeVar\howto</a>              - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\HighDimensionalModelRepresentation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\HighDimensionalModelRepresentation\howto">examples\HighDimensionalModelRepresentation\howto</a>                  - % housekeeping
%   <a href="matlab:help examples\HighDimensionalModelRepresentation\ishigami">examples\HighDimensionalModelRepresentation\ishigami</a>               - %%========================================
%   <a href="matlab:help examples\HighDimensionalModelRepresentation\polynomial_integration">examples\HighDimensionalModelRepresentation\polynomial_integration</a> -  is of the form a0+a1*x+...+ar*x^r
%   examples\HighDimensionalModelRepresentation\satelli_sobol95        - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\HybridExpectations   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\HybridExpectations\howto">examples\HybridExpectations\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\LagsOnExogenousProcesses   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\LagsOnExogenousProcesses\howto">examples\LagsOnExogenousProcesses\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\Linear_vs_Loglinear   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\Linear_vs_Loglinear\howto">examples\Linear_vs_Loglinear\howto</a>                           - % housekeeping
%   <a href="matlab:help examples\Linear_vs_Loglinear\tshocksnk_steadystate">examples\Linear_vs_Loglinear\tshocksnk_steadystate</a>           - params=vertcat(param_obj.startval); %#ok<NASGU>
%   <a href="matlab:help examples\Linear_vs_Loglinear\tshocksnk_steadystate_loglinear">examples\Linear_vs_Loglinear\tshocksnk_steadystate_loglinear</a> - params=vertcat(param_obj.startval); %#ok<NASGU>
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\BiTraum2013   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\Untitled5">examples\MarkovSwitching\BiTraum2013\Untitled5</a>               - percent sign','%'
%   examples\MarkovSwitching\BiTraum2013\all_items               - (No help available)
%   examples\MarkovSwitching\BiTraum2013\bi_traum_13_ssfile      - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\bt_ssfile">examples\MarkovSwitching\BiTraum2013\bt_ssfile</a>               - Instruct RISE not to bother checking that this is the true steady state
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\dynamic_derivatives">examples\MarkovSwitching\BiTraum2013\dynamic_derivatives</a>     - anonymous function for checking validity
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\endo_exo">examples\MarkovSwitching\BiTraum2013\endo_exo</a>                - % Code automagically generated by rise on 13-Jul-2013 13:47:20
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\finding_bounds_for_rise">examples\MarkovSwitching\BiTraum2013\finding_bounds_for_rise</a> - % housekeeping
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\first_order_solver">examples\MarkovSwitching\BiTraum2013\first_order_solver</a>      -  options are 0 or 'msre_klein' (for constant-parameter models)
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\junk_print">examples\MarkovSwitching\BiTraum2013\junk_print</a>              - prototype{3}=vi.obj.ncalls;
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\msre_filter_kim_nelson">examples\MarkovSwitching\BiTraum2013\msre_filter_kim_nelson</a>  - this filter assumes a state space of the form
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\my_print">examples\MarkovSwitching\BiTraum2013\my_print</a>                - min_ncount) % minimum number of occurrences at which substitution occurs
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\myderivs">examples\MarkovSwitching\BiTraum2013\myderivs</a>                - % Code automagically generated by rise on 27-Jul-2013 21:50:07
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\parse_declaration">examples\MarkovSwitching\BiTraum2013\parse_declaration</a>       - classdef parse_declaration < handle
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\planner_stuff">examples\MarkovSwitching\BiTraum2013\planner_stuff</a>           - % Code automagically generated by rise on 14-Jul-2013 07:13:15
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\simulate_second_order">examples\MarkovSwitching\BiTraum2013\simulate_second_order</a>   - Q=solution.Q;
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\state_matrices">examples\MarkovSwitching\BiTraum2013\state_matrices</a>          - examples of calls:
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\taoMacProblem">examples\MarkovSwitching\BiTraum2013\taoMacProblem</a>           - housekeeping
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\transit_matrix">examples\MarkovSwitching\BiTraum2013\transit_matrix</a>          - % Code automagically generated by rise on 05-Jul-2013 10:50:41
%   <a href="matlab:help examples\MarkovSwitching\BiTraum2013\transition_matrix">examples\MarkovSwitching\BiTraum2013\transition_matrix</a>       - % Code automagically generated by rise on 21-Jul-2013 16:12:48
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\FarmerWaggonerZha2010   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\FarmerWaggonerZha2010\fwz_test_suite_dynare_conf">examples\MarkovSwitching\FarmerWaggonerZha2010\fwz_test_suite_dynare_conf</a>                              - % housekeeping
%   examples\MarkovSwitching\FarmerWaggonerZha2010\generate_markov_switching_rational_expectations_problem - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\FarmerWaggonerZha2010\howto">examples\MarkovSwitching\FarmerWaggonerZha2010\howto</a>                                                   - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\GavinKeenRichterThrockmorton   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant</a>                         - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_dynamic">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_dynamic</a>                 - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_set_auxiliary_variables">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_set_auxiliary_variables</a> - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_static">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_static</a>                  - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_steadystate2">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant_steadystate2</a>            -  state generated by Dynare preprocessor
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\howto">examples\MarkovSwitching\GavinKeenRichterThrockmorton\howto</a>                            - % housekeeping
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\planar_diff">examples\MarkovSwitching\GavinKeenRichterThrockmorton\planar_diff</a>                      - Current=cell(2,incmnt); % 1- indexes, 2-pointers
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\stoch_simul____">examples\MarkovSwitching\GavinKeenRichterThrockmorton\stoch_simul____</a>                  - should also allow for passing options
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Endogenous_Shocks_Parameters___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Endogenous_Shocks_Parameters___</a> - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Parameters___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Parameters___</a>                   - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\StaticEndogenous___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\StaticEndogenous___</a>             - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Static_BGP_Endogenous___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\Static_BGP_Endogenous___</a>        - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\balanced_growth___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\balanced_growth___</a>              - % Code automagically generated by rise on 04-Sep-2013 09:46:37
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\definitions___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\definitions___</a>                  - % Code automagically generated by rise on 04-Sep-2013 09:46:37
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\dynamic___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\dynamic___</a>                      - % Code automagically generated by rise on 04-Sep-2013 09:46:37
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\dynamic_params___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\dynamic_params___</a>               - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\static___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\static___</a>                       - % Code automagically generated by rise on 04-Sep-2013 09:46:37
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\steady_state_model___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\steady_state_model___</a>           - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\transition_matrix___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\transition_matrix___</a>            - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\vectorized_dynamic___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\vectorized_dynamic___</a>           - % Code automagically generated by rise on 04-Sep-2013 09:46:37
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\vectorized_dynamic_params___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\constant\routines\vectorized_dynamic_params___</a>    - % Code automagically generated by rise on 04-Sep-2013 09:46:38
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Endogenous_Shocks_Parameters___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Endogenous_Shocks_Parameters___</a> - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Parameters___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Parameters___</a>                   - % Code automagically generated by rise on 02-Sep-2013 14:51:53
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\StaticEndogenous___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\StaticEndogenous___</a>             - % Code automagically generated by rise on 02-Sep-2013 14:51:53
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Static_BGP_Endogenous___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\Static_BGP_Endogenous___</a>        - % Code automagically generated by rise on 02-Sep-2013 14:51:53
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\balanced_growth___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\balanced_growth___</a>              - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\definitions___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\definitions___</a>                  - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\dynamic___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\dynamic___</a>                      - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\dynamic_params___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\dynamic_params___</a>               - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\static___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\static___</a>                       - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\steady_state_model___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\steady_state_model___</a>           - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\transition_matrix___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\transition_matrix___</a>            - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\vectorized_dynamic___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\vectorized_dynamic___</a>           - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%   <a href="matlab:help examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\vectorized_dynamic_params___">examples\MarkovSwitching\GavinKeenRichterThrockmorton\globdynzlb\routines\vectorized_dynamic_params___</a>    - % Code automagically generated by rise on 02-Sep-2013 14:51:52
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial1   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial1\howto">examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial1\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial2   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial2\howto">examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial2\howto</a>        - %
%   <a href="matlab:help examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial2\howto_report">examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial2\howto_report</a> - housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial3   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial3\howto">examples\MarkovSwitching\LiuWaggonerZha2009\Tutorial3\howto</a> - %
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\NorwayMonetaryPolicy   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\NorwayMonetaryPolicy\RegimeSwitchesNorway">examples\MarkovSwitching\NorwayMonetaryPolicy\RegimeSwitchesNorway</a>           - % housekeeping
%   <a href="matlab:help examples\MarkovSwitching\NorwayMonetaryPolicy\howto_forecast">examples\MarkovSwitching\NorwayMonetaryPolicy\howto_forecast</a>                 - % introduction
%   <a href="matlab:help examples\MarkovSwitching\NorwayMonetaryPolicy\manhattan_beach">examples\MarkovSwitching\NorwayMonetaryPolicy\manhattan_beach</a>                - % housekeeping
%   <a href="matlab:help examples\MarkovSwitching\NorwayMonetaryPolicy\max_operator_behavior_test">examples\MarkovSwitching\NorwayMonetaryPolicy\max_operator_behavior_test</a>     - this function checks that expected monetary policy shocks are positive
%   <a href="matlab:help examples\MarkovSwitching\NorwayMonetaryPolicy\steady_state_4_Canonical_Const">examples\MarkovSwitching\NorwayMonetaryPolicy\steady_state_4_Canonical_Const</a> - unction [ss,var_names,retcode]=steady_state_4_Canonical_Const(param_obj)%param_struct=struct(parlist,parvals)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\RamirezWaggonerZha2012   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare</a>                              - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_dynamic">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_dynamic</a>                      - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_set_auxiliary_variables">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_set_auxiliary_variables</a>      - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_static">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_static</a>                       - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_steadystate2">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_steadystate2</a>                 -  state generated by Dynare preprocessor
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_test_derivs">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_test_derivs</a>                         - %
%   examples\MarkovSwitching\RamirezWaggonerZha2012\functional_iteration_convergence_conditions - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\howto">examples\MarkovSwitching\RamirezWaggonerZha2012\howto</a>                                       - % Housekeeping
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\howto_frwz_nk">examples\MarkovSwitching\RamirezWaggonerZha2012\howto_frwz_nk</a>                               - % Housekeeping
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\rwz_steady_state">examples\MarkovSwitching\RamirezWaggonerZha2012\rwz_steady_state</a>                            - unction [params,ss,retcode,imposed]=rwz_steady_state(params,flag)%param_struct=struct(parlist,parvals)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Endogenous_Shocks_Parameters___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Endogenous_Shocks_Parameters___</a> - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Parameters___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Parameters___</a>                   - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\StaticEndogenous___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\StaticEndogenous___</a>             - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Static_BGP_Endogenous___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\Static_BGP_Endogenous___</a>        - % Code automagically generated by rise on 20-Aug-2013 15:52:37
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\balanced_growth___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\balanced_growth___</a>              - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\definitions___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\definitions___</a>                  - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\dynamic___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\dynamic___</a>                      - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\dynamic_params___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\dynamic_params___</a>               - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\static___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\static___</a>                       - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\steady_state_model___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\steady_state_model___</a>           - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\transition_matrix___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\transition_matrix___</a>            - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\vectorized_dynamic___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\vectorized_dynamic___</a>           - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\vectorized_dynamic_params___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk\routines\vectorized_dynamic_params___</a>    - % Code automagically generated by rise on 20-Aug-2013 15:52:36
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Endogenous_Shocks_Parameters___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Endogenous_Shocks_Parameters___</a> - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Parameters___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Parameters___</a>                   - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\StaticEndogenous___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\StaticEndogenous___</a>             - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Static_BGP_Endogenous___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\Static_BGP_Endogenous___</a>        - % Code automagically generated by rise on 14-Aug-2013 12:04:44
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\balanced_growth___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\balanced_growth___</a>              - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\definitions___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\definitions___</a>                  - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\dynamic___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\dynamic___</a>                      - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\dynamic_params___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\dynamic_params___</a>               - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\static___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\static___</a>                       - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\steady_state_model___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\steady_state_model___</a>           - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\transition_matrix___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\transition_matrix___</a>            - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\vectorized_dynamic___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\vectorized_dynamic___</a>           - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\vectorized_dynamic_params___">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_nk_dynare_test\routines\vectorized_dynamic_params___</a>    - % Code automagically generated by rise on 14-Aug-2013 12:04:43
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\ModelWithMeasurementErrors   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\ModelWithMeasurementErrors\howto">examples\ModelWithMeasurementErrors\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\ModelsWithSteadyStateFile   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\ModelsWithSteadyStateFile\howto">examples\ModelsWithSteadyStateFile\howto</a>                          - % housekeeping
%   <a href="matlab:help examples\ModelsWithSteadyStateFile\steady_state_4_Canonical_Const">examples\ModelsWithSteadyStateFile\steady_state_4_Canonical_Const</a> - unction [ss,param_obj,retcode,imposed]=steady_state_4_Canonical_Const(param_obj,flag)%param_struct=struct(parlist,parvals)
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MomentsToBounds   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MomentsToBounds\howto">examples\MomentsToBounds\howto</a> - % computing the bounds from the mean and standard deviations from dynare
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MonteCarloFiltering   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MonteCarloFiltering\howto">examples\MonteCarloFiltering\howto</a>        - % housekeeping
%   <a href="matlab:help examples\MonteCarloFiltering\price_puzzle">examples\MonteCarloFiltering\price_puzzle</a> - % inflation
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\MonteCarloFiltering\trinity\routines   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\Endogenous_Shocks_Parameters___">examples\MonteCarloFiltering\trinity\routines\Endogenous_Shocks_Parameters___</a> - % Code automagically generated by rise on 06-Sep-2013 13:04:56
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\Parameters___">examples\MonteCarloFiltering\trinity\routines\Parameters___</a>                   - % Code automagically generated by rise on 06-Sep-2013 13:04:56
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\StaticEndogenous___">examples\MonteCarloFiltering\trinity\routines\StaticEndogenous___</a>             - % Code automagically generated by rise on 06-Sep-2013 13:04:56
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\Static_BGP_Endogenous___">examples\MonteCarloFiltering\trinity\routines\Static_BGP_Endogenous___</a>        - % Code automagically generated by rise on 06-Sep-2013 13:04:56
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\balanced_growth___">examples\MonteCarloFiltering\trinity\routines\balanced_growth___</a>              - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\definitions___">examples\MonteCarloFiltering\trinity\routines\definitions___</a>                  - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\dynamic___">examples\MonteCarloFiltering\trinity\routines\dynamic___</a>                      - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\dynamic_params___">examples\MonteCarloFiltering\trinity\routines\dynamic_params___</a>               - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\static___">examples\MonteCarloFiltering\trinity\routines\static___</a>                       - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\steady_state_model___">examples\MonteCarloFiltering\trinity\routines\steady_state_model___</a>           - % Code automagically generated by rise on 06-Sep-2013 13:04:56
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\transition_matrix___">examples\MonteCarloFiltering\trinity\routines\transition_matrix___</a>            - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\vectorized_dynamic___">examples\MonteCarloFiltering\trinity\routines\vectorized_dynamic___</a>           - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%   <a href="matlab:help examples\MonteCarloFiltering\trinity\routines\vectorized_dynamic_params___">examples\MonteCarloFiltering\trinity\routines\vectorized_dynamic_params___</a>    - % Code automagically generated by rise on 06-Sep-2013 13:04:55
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\NonlinearModels   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\NonlinearModels\fs2000_steadystate">examples\NonlinearModels\fs2000_steadystate</a>         - computes the steady state of fs2000 analyticaly
%   <a href="matlab:help examples\NonlinearModels\fs2000_steadystate_initval">examples\NonlinearModels\fs2000_steadystate_initval</a> - gives approximate values for the steady state of fs2000
%   <a href="matlab:help examples\NonlinearModels\howto">examples\NonlinearModels\howto</a>                      - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\OptimalSimpleRules   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\OptimalSimpleRules\howto">examples\OptimalSimpleRules\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\SVAR   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\SVAR\SvarTutorial">examples\SVAR\SvarTutorial</a>                   - % housekeeping
%   <a href="matlab:help examples\SVAR\blanchard_perotti_2002">examples\SVAR\blanchard_perotti_2002</a>         - % housekeeping
%   <a href="matlab:help examples\SVAR\blanchard_perotti_restrictions">examples\SVAR\blanchard_perotti_restrictions</a> - a_g=param(1); % -inf inf
%   <a href="matlab:help examples\SVAR\mertens_ravn_data">examples\SVAR\mertens_ravn_data</a>              -  used in papers:
%   <a href="matlab:help examples\SVAR\peersman_data">examples\SVAR\peersman_data</a>                  - source: http://www.cambridge.org/features/econmodelling/download/MATLAB.zip
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\StickyInformation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\StickyInformation\howto">examples\StickyInformation\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\StochasticReplanning_switching_Estimation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\howto">examples\StochasticReplanning_switching_Estimation\howto</a>                   - % housekeeping
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\howto_steady_state">examples\StochasticReplanning_switching_Estimation\howto_steady_state</a>      - % housekeeping
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\usmodel_lc_steady_state">examples\StochasticReplanning_switching_Estimation\usmodel_lc_steady_state</a> - % push back the parameter inside the params
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\StochasticReplanning_switching_Estimation\archive   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\archive\usmodel_steadystate">examples\StochasticReplanning_switching_Estimation\archive\usmodel_steadystate</a> - computes the steady state for the observed variables in the smets-wouters
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\SymbolicDifferentiation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\SymbolicDifferentiation\test_symbolic">examples\SymbolicDifferentiation\test_symbolic</a>   - % housekeeping
%   <a href="matlab:help examples\SymbolicDifferentiation\test_symbolic_2">examples\SymbolicDifferentiation\test_symbolic_2</a> - % create function
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\TimeSeriesObjects   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\TimeSeriesObjects\howto">examples\TimeSeriesObjects\howto</a> - % create empty time series
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\KaijiTutorial   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\VariousModels\KaijiTutorial\howto">examples\VariousModels\KaijiTutorial\howto</a>     - % housekeeping
%   <a href="matlab:help examples\VariousModels\KaijiTutorial\howto_exp">examples\VariousModels\KaijiTutorial\howto_exp</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\LiuWangZha2009   %%%%%%%%%%%%%%%%%%%%
%
%   examples\VariousModels\LiuWangZha2009\data_diff_FHFA - (No help available)
%   <a href="matlab:help examples\VariousModels\LiuWangZha2009\howto">examples\VariousModels\LiuWangZha2009\howto</a>          - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\PeterIreland\technology_shocks   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\VariousModels\PeterIreland\technology_shocks\howto">examples\VariousModels\PeterIreland\technology_shocks\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\SmetsWouters   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\VariousModels\SmetsWouters\howto">examples\VariousModels\SmetsWouters\howto</a>               - % housekeeping
%   <a href="matlab:help examples\VariousModels\SmetsWouters\usmodel_steadystate">examples\VariousModels\SmetsWouters\usmodel_steadystate</a> - computes the steady state for the observed variables in the smets-wouters
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\fs2000   %%%%%%%%%%%%%%%%%%%%
%
%   examples\VariousModels\fs2000\csminwellwrap - (No help available)
%   <a href="matlab:help examples\VariousModels\fs2000\howto">examples\VariousModels\fs2000\howto</a>         - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\VariousModels\fs2000\dynare_version   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\VariousModels\fs2000\dynare_version\dsge_likelihood">examples\VariousModels\fs2000\dynare_version\dsge_likelihood</a>       - Evaluates the posterior kernel of a dsge model.
%   <a href="matlab:help examples\VariousModels\fs2000\dynare_version\evaluate_steady_state">examples\VariousModels\fs2000\dynare_version\evaluate_steady_state</a> - function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
%   <a href="matlab:help examples\VariousModels\fs2000\dynare_version\fsdat_simul">examples\VariousModels\fs2000\dynare_version\fsdat_simul</a>           - Generated data, used by fs2000.mod
%   <a href="matlab:help examples\VariousModels\fs2000\dynare_version\priordens">examples\VariousModels\fs2000\dynare_version\priordens</a>             - Computes a prior density for the structural parameters of DSGE models
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\differentiation   %%%%%%%%%%%%%%%%%%%%
%
%   m\differentiation\rise_numjac - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\filtering   %%%%%%%%%%%%%%%%%%%%
%
%   m\filtering\conditional_likelihood           - (No help available)
%   m\filtering\initial_markov_distribution      - (No help available)
%   <a href="matlab:help m\filtering\kalman_initialization">m\filtering\kalman_initialization</a>            - There is no documentation of this function yet.
%   <a href="matlab:help m\filtering\likelihood_dsge_var">m\filtering\likelihood_dsge_var</a>              - N.B: This function assumes the likelihood is to be maximized
%   <a href="matlab:help m\filtering\likelihood_markov_switching_dsge">m\filtering\likelihood_markov_switching_dsge</a> - unction [LogLik,Incr,retcode,obj]=likelihood_markov_switching_dsge(params,obj)%
%   <a href="matlab:help m\filtering\likelihood_optimal_simple_rule">m\filtering\likelihood_optimal_simple_rule</a>   - % this important output is not created yet
%   <a href="matlab:help m\filtering\msre_kalman_cell">m\filtering\msre_kalman_cell</a>                 - this filter assumes a state space of the form
%   <a href="matlab:help m\filtering\msre_linear_filter">m\filtering\msre_linear_filter</a>               - all rows of Q should sum to 1
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\filtering\+junk   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\filtering\+junk\CheckPositiveDefiniteness">m\filtering\+junk\CheckPositiveDefiniteness</a>                - iF=F\eye(size(F,1));%inv(F);
%   <a href="matlab:help m\filtering\+junk\kalman_prediction">m\filtering\+junk\kalman_prediction</a>                        - only compute the places where there is some action
%   <a href="matlab:help m\filtering\+junk\kalman_update">m\filtering\+junk\kalman_update</a>                            - K=P(:,obs_id);
%   <a href="matlab:help m\filtering\+junk\markov_switching_kalman_filter">m\filtering\+junk\markov_switching_kalman_filter</a>           - Detailed explanation to come here
%   <a href="matlab:help m\filtering\+junk\markov_switching_kalman_filter_real_time">m\filtering\+junk\markov_switching_kalman_filter_real_time</a> - y,... % data
%   <a href="matlab:help m\filtering\+junk\smoothing_step">m\filtering\+junk\smoothing_step</a>                           - Note that Durbin and Koopman define K=T*P*Z'*iF, while here it is defined
%   <a href="matlab:help m\filtering\+junk\symmetrize">m\filtering\+junk\symmetrize</a>                               - unction B=symmetrize(A)%,flag
%   <a href="matlab:help m\filtering\+junk\update_and_collapse">m\filtering\+junk\update_and_collapse</a>                      - Ptt(:,:,snow)=Ptt(:,:,snow)/PAItt(snow); %symmetrize()
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\forecasting   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\forecasting\forecast_engine">m\forecasting\forecast_engine</a>              - Description:
%   <a href="matlab:help m\forecasting\forecasting_engine">m\forecasting\forecasting_engine</a>           - Description:
%   <a href="matlab:help m\forecasting\three_pass_regression_filter">m\forecasting\three_pass_regression_filter</a> - inputs:
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\forecasting\forecasting_tools   %%%%%%%%%%%%%%%%%%%%
%
%   m\forecasting\forecasting_tools\BuildShockRestrictions                         - (No help available)
%   <a href="matlab:help m\forecasting\forecasting_tools\ComputeForecasts">m\forecasting\forecasting_tools\ComputeForecasts</a>                               -  conditional forecasts
%   <a href="matlab:help m\forecasting\forecasting_tools\ConditionalDistributionOfEndogenousAndShocks">m\forecasting\forecasting_tools\ConditionalDistributionOfEndogenousAndShocks</a>   - =============================
%   <a href="matlab:help m\forecasting\forecasting_tools\ConditionalProjectionSubEngine">m\forecasting\forecasting_tools\ConditionalProjectionSubEngine</a>                 - EndogenousConditions,ShocksConditions,verbose) % ,Yf,
%   <a href="matlab:help m\forecasting\forecasting_tools\ConditionalStateMatrices">m\forecasting\forecasting_tools\ConditionalStateMatrices</a>                       - indicate that an element is a covariance matrix by putting it into a cell
%   m\forecasting\forecasting_tools\InitializeArray                                - (No help available)
%   <a href="matlab:help m\forecasting\forecasting_tools\NullAndColumnSpaces">m\forecasting\forecasting_tools\NullAndColumnSpaces</a>                            - 2=V(:,1:q); % another basis which is different but gives the same end results is M22=null(M1')
%   <a href="matlab:help m\forecasting\forecasting_tools\RemoveHoles">m\forecasting\forecasting_tools\RemoveHoles</a>                                    - % check whether V is a covariance matrix
%   <a href="matlab:help m\forecasting\forecasting_tools\TruncatedMultivariateNormalRnd">m\forecasting\forecasting_tools\TruncatedMultivariateNormalRnd</a>                 - x=TruncatedMultivariateNormalRnd(mu,SIG,lb,ub)
%   m\forecasting\forecasting_tools\UnconditionalDistributionOfEndogenousAndShocks - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\hessian_computation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\hessian_computation\finite_difference_hessian">m\hessian_computation\finite_difference_hessian</a> - Compute the stepsize (h)
%   <a href="matlab:help m\hessian_computation\hessian_computation_test">m\hessian_computation\hessian_computation_test</a>  - load some data
%   <a href="matlab:help m\hessian_computation\hessian_conditioner">m\hessian_computation\hessian_conditioner</a>       - based on ddom.m
%   <a href="matlab:help m\hessian_computation\outer_product_hessian">m\hessian_computation\outer_product_hessian</a>     - Objective is assumed to have at least two output arguments
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\optimizers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\optimizers\bat">m\optimizers\bat</a>                          - Reference:
%   <a href="matlab:help m\optimizers\batm">m\optimizers\batm</a>                         - Reference:
%   <a href="matlab:help m\optimizers\bbo_3">m\optimizers\bbo_3</a>                        - 'mutation_type','gaussian' % 'cauchy','uniform'
%   <a href="matlab:help m\optimizers\bbo_gate">m\optimizers\bbo_gate</a>                     - Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   <a href="matlab:help m\optimizers\bbo_gate_conclude">m\optimizers\bbo_gate_conclude</a>            - exitflag=0; % disp('Too many function evaluations or iterations.')
%   <a href="matlab:help m\optimizers\bee">m\optimizers\bee</a>                          - lassdef bee %< handle
%   <a href="matlab:help m\optimizers\bee_2">m\optimizers\bee_2</a>                        - Reference: Inspired from Karaboga
%   <a href="matlab:help m\optimizers\bee__">m\optimizers\bee__</a>                        - this function sorts using f while the class alternative sorts using
%   <a href="matlab:help m\optimizers\bee_demo">m\optimizers\bee_demo</a>                     - % optimizer-specific properties
%   <a href="matlab:help m\optimizers\bee_demo_launch">m\optimizers\bee_demo_launch</a>              - =2; % number of parameters
%   <a href="matlab:help m\optimizers\bee_gate">m\optimizers\bee_gate</a>                     -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\bee_gate_conclude">m\optimizers\bee_gate_conclude</a>            - exitflag=0; % disp('Too many function evaluations or iterations.')
%   <a href="matlab:help m\optimizers\blockwise_optimization">m\optimizers\blockwise_optimization</a>       - get the name of the optimizer and check whether it comes from matlab, in
%   <a href="matlab:help m\optimizers\cmsa">m\optimizers\cmsa</a>                         - % algorithm specific options
%   <a href="matlab:help m\optimizers\cmsa_gate">m\optimizers\cmsa_gate</a>                    -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\cmsa_gate_conclude">m\optimizers\cmsa_gate_conclude</a>           - exitflag=0; % disp('Too many function evaluations or iterations.')
%   m\optimizers\csminwellwrap                - (No help available)
%   <a href="matlab:help m\optimizers\diff_lev_bat">m\optimizers\diff_lev_bat</a>                 - Reference:
%   <a href="matlab:help m\optimizers\evaluate_individual">m\optimizers\evaluate_individual</a>          - % correct the bounds
%   <a href="matlab:help m\optimizers\gampc">m\optimizers\gampc</a>                        - % algorithm specific options
%   <a href="matlab:help m\optimizers\gampc_2">m\optimizers\gampc_2</a>                      -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\hgpsal">m\optimizers\hgpsal</a>                       - If just 'defaults' passed in, return the default options in X
%   <a href="matlab:help m\optimizers\hybrid_artificial_bee_colony">m\optimizers\hybrid_artificial_bee_colony</a> - 'F',0.5 % F=2*rand; %[0,2]
%   <a href="matlab:help m\optimizers\iwo">m\optimizers\iwo</a>                          - IWO: Invasive Weed Optimization
%   <a href="matlab:help m\optimizers\local_optimize">m\optimizers\local_optimize</a>               - trim everything
%   <a href="matlab:help m\optimizers\local_optimize_2">m\optimizers\local_optimize_2</a>             - trim everything
%   <a href="matlab:help m\optimizers\local_optimize_gate">m\optimizers\local_optimize_gate</a>          - if obj.iterations>=obj.MaxIter || ...
%   <a href="matlab:help m\optimizers\local_reoptimize">m\optimizers\local_reoptimize</a>             - trim everything
%   <a href="matlab:help m\optimizers\mbo">m\optimizers\mbo</a>                          - migrating birds optimization (MBO)
%   <a href="matlab:help m\optimizers\msnlp">m\optimizers\msnlp</a>                        - If just 'defaults' passed in, return the default options in X
%   <a href="matlab:help m\optimizers\optimtestfun">m\optimizers\optimtestfun</a>                 - function u=u_func(x,a,k,m)
%   m\optimizers\rebuild_population           - (No help available)
%   <a href="matlab:help m\optimizers\rise_nelder_mead">m\optimizers\rise_nelder_mead</a>             - simplex algorithm with nonlinear constraints
%   <a href="matlab:help m\optimizers\studga">m\optimizers\studga</a>                       - % optimizer-specific options
%   <a href="matlab:help m\optimizers\studga_gate">m\optimizers\studga_gate</a>                  -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\studga_gate_conclude">m\optimizers\studga_gate_conclude</a>         - exitflag=0; % disp('Too many function evaluations or iterations.')
%   m\optimizers\studga_gate_iris             - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\optimizers\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\optimizers\private\check_convergence">m\optimizers\private\check_convergence</a>              - if rem(obj.iterations,100)==0 || ~isempty(stopflag)
%   <a href="matlab:help m\optimizers\private\clear_duplicates">m\optimizers\private\clear_duplicates</a>               - % select the chromosomes to change randomly
%   m\optimizers\private\compare_individuals            - (No help available)
%   m\optimizers\private\compute_fitness                - (No help available)
%   m\optimizers\private\dispersion                     - (No help available)
%   <a href="matlab:help m\optimizers\private\display_progress">m\optimizers\private\display_progress</a>               - fprintf(1,'restart # %3.0f   iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    stdev %8.4f    f-Count  %8.0f   routine %s\n',...
%   m\optimizers\private\distance                       - (No help available)
%   m\optimizers\private\dynamic_penalty                - (No help available)
%   m\optimizers\private\find_farthest                  - (No help available)
%   m\optimizers\private\find_nearest                   - (No help available)
%   <a href="matlab:help m\optimizers\private\generate_candidates">m\optimizers\private\generate_candidates</a>            -  we get 4 outputs?
%   <a href="matlab:help m\optimizers\private\manual_stopping">m\optimizers\private\manual_stopping</a>                - rawfile = char(textread(ManualStoppingFile,'%s','delimiter','\n','whitespace','','bufsize',40000));
%   <a href="matlab:help m\optimizers\private\optimization_universal_options">m\optimizers\private\optimization_universal_options</a> - optimization_universal_options: sets the common default options for all
%   m\optimizers\private\recenter                       - (No help available)
%   m\optimizers\private\selection_process              - (No help available)
%   <a href="matlab:help m\optimizers\private\sort_population">m\optimizers\private\sort_population</a>                - % Deb sorting
%   <a href="matlab:help m\optimizers\private\uniform_sampling">m\optimizers\private\uniform_sampling</a>               - samples without repetition
%   <a href="matlab:help m\optimizers\private\weighted_sampling">m\optimizers\private\weighted_sampling</a>              - samples without repetition
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\parsers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\parsers\analytical_symbolic_form">m\parsers\analytical_symbolic_form</a>       - turns x(1),y(1),ss(1),p(1) into x_1,y_1,ss_1,p_1
%   <a href="matlab:help m\parsers\cellstr2mat">m\parsers\cellstr2mat</a>                    - unction out=cellstr2mat(CellMat) %,MatName
%   <a href="matlab:help m\parsers\chain_grid">m\parsers\chain_grid</a>                     - v is a vector of the number of states in each markov chain
%   <a href="matlab:help m\parsers\collect_symbolic_list">m\parsers\collect_symbolic_list</a>          - list_ss_=[list_ss_,list_ss{ilist}]; %#ok<AGROW>
%   m\parsers\find_occurrences               - (No help available)
%   <a href="matlab:help m\parsers\find_separators">m\parsers\find_separators</a>                - location of delimiters
%   m\parsers\is_atom                        - (No help available)
%   <a href="matlab:help m\parsers\match_parentheses">m\parsers\match_parentheses</a>              - if left_right(ii)==left_par %<--- strcmp(left_right(ii),left_par)
%   <a href="matlab:help m\parsers\my_strtok">m\parsers\my_strtok</a>                      - delims=eqtn(1:start-1)
%   <a href="matlab:help m\parsers\parenthesize">m\parsers\parenthesize</a>                   - func is any of the following: '^', '/' , '-', '+' ,'*'
%   <a href="matlab:help m\parsers\remove_unnecessary_parentheses">m\parsers\remove_unnecessary_parentheses</a> -  parentheses around atoms and multiplicative expressions
%   <a href="matlab:help m\parsers\rise_algebra_cat">m\parsers\rise_algebra_cat</a>               - % do nothing
%   <a href="matlab:help m\parsers\rise_isa">m\parsers\rise_isa</a>                       - rise_isa(string) determines the types of the string
%   <a href="matlab:help m\parsers\symbolic2model">m\parsers\symbolic2model</a>                 - shad=[shad,eqtn]; %#ok<*AGROW>
%   <a href="matlab:help m\parsers\trim_symbolic_equation">m\parsers\trim_symbolic_equation</a>         - %
%   <a href="matlab:help m\parsers\update_markov_chains_info">m\parsers\update_markov_chains_info</a>      - creates a markov chain info and updates it as new information comes in.
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\parsers\+string_optimize   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\parsers\+string_optimize\remove_unnecessary_parentheses">m\parsers\+string_optimize\remove_unnecessary_parentheses</a> - Replace "((....))" with "(....)"
%   <a href="matlab:help m\parsers\+string_optimize\remove_zeros_and_ones">m\parsers\+string_optimize\remove_zeros_and_ones</a>          - % zeros in loops
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\solvers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\solvers\collapse_array">m\solvers\collapse_array</a>                               - takes as input a multidimensional array and outputs a cell array
%   <a href="matlab:help m\solvers\computational_savings">m\solvers\computational_savings</a>                        - this function separates static variables from dynamic ones. It places all
%   <a href="matlab:help m\solvers\dsge_lc_solve">m\solvers\dsge_lc_solve</a>                                - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\dsge_solve_aim">m\solvers\dsge_solve_aim</a>                               - ags=1; % no of lags and leads
%   <a href="matlab:help m\solvers\dsge_solve_gensys">m\solvers\dsge_solve_gensys</a>                            - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\dsge_solve_klein">m\solvers\dsge_solve_klein</a>                             - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\expand_array">m\solvers\expand_array</a>                                 - case 1	% A0, Aminus or Aplus
%   <a href="matlab:help m\solvers\fix_point_iterator">m\solvers\fix_point_iterator</a>                           - this function solves for a fix point. Inputs are:
%   <a href="matlab:help m\solvers\gensys">m\solvers\gensys</a>                                       - function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
%   <a href="matlab:help m\solvers\get_default_optimization_option">m\solvers\get_default_optimization_option</a>              - lc_reconvexify:
%   m\solvers\is_eigenvalue_solver_candidate               - (No help available)
%   <a href="matlab:help m\solvers\loose_commitment_solver">m\solvers\loose_commitment_solver</a>                      - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\loose_commitment_solver_fix_point_unfinished">m\solvers\loose_commitment_solver_fix_point_unfinished</a> - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\loose_commitment_to_markov_switching">m\solvers\loose_commitment_to_markov_switching</a>         - this function puts the loose commitment solution into a markov switching
%   m\solvers\markov_switching_dsge_objective              - (No help available)
%   <a href="matlab:help m\solvers\markov_switching_dsge_stack">m\solvers\markov_switching_dsge_stack</a>                  - C_st,... % endo_nbr x h matrix of constant
%   <a href="matlab:help m\solvers\msre_aim">m\solvers\msre_aim</a>                                     - solve for T only
%   <a href="matlab:help m\solvers\msre_gensys">m\solvers\msre_gensys</a>                                  - solve for T only
%   <a href="matlab:help m\solvers\msre_klein">m\solvers\msre_klein</a>                                   - solve for T only
%   <a href="matlab:help m\solvers\msre_matrix_times_vector">m\solvers\msre_matrix_times_vector</a>                     - the old version is faster and solves but solves the problem
%   <a href="matlab:help m\solvers\msre_solve">m\solvers\msre_solve</a>                                   - This procedure assumes the steady state has been solved and that apart
%   <a href="matlab:help m\solvers\qzdiv">m\solvers\qzdiv</a>                                        - function [A,B,Q,Z] = qzdiv(stake,A,B,Q,Z)
%   <a href="matlab:help m\solvers\qzswitch">m\solvers\qzswitch</a>                                     - function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%   <a href="matlab:help m\solvers\schur_solver">m\solvers\schur_solver</a>                                 - T,S,Q,Z] = qz(F,G,'real');%complex
%   <a href="matlab:help m\solvers\solve_steady_state">m\solvers\solve_steady_state</a>                           - %         ys0=0*ys0;
%   <a href="matlab:help m\solvers\transpose_free_quasi_minimum_residual">m\solvers\transpose_free_quasi_minimum_residual</a>        - A,... % coefficient matrix
%   <a href="matlab:help m\solvers\unearth_frwrd_matrix">m\solvers\unearth_frwrd_matrix</a>                         - this function extracts Aplus from Gplus and is used for solving models
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\solvers\X_equal_A_X_B_plus_C_solvers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\diffuse_lyapunov_equation">m\solvers\X_equal_A_X_B_plus_C_solvers\diffuse_lyapunov_equation</a> - % check whether there are diffuse_all elements as they require special treatment
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\doubling_solve">m\solvers\X_equal_A_X_B_plus_C_solvers\doubling_solve</a>            - unction [P,retcode]=doubling_solve(A,B,C,options)%MaxIter,TolFun,verbose
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\lyapunov_equation">m\solvers\X_equal_A_X_B_plus_C_solvers\lyapunov_equation</a>         - example
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_a_la_tadonki">m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_a_la_tadonki</a>     - attempts to solve the equation V=A*V*B+C
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_solve">m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_solve</a>            - solves the linear equation X=A*X*B+C
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\tests">m\solvers\X_equal_A_X_B_plus_C_solvers\tests</a>                     - %
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\utilities   %%%%%%%%%%%%%%%%%%%%
%
%   m\utilities\CheckArgument                            - (No help available)
%   m\utilities\alpha_probability                        - (No help available)
%   <a href="matlab:help m\utilities\build_grid">m\utilities\build_grid</a>                               - see also mygrid
%   <a href="matlab:help m\utilities\cell2object">m\utilities\cell2object</a>                              - ypes=C(3,:); %#ok<NASGU>
%   <a href="matlab:help m\utilities\code2file">m\utilities\code2file</a>                                - writes code in char form to an m-file, which could be a function if
%   <a href="matlab:help m\utilities\concatenate_series_from_different_models">m\utilities\concatenate_series_from_different_models</a> - dbcell is a cell array of structures of time series objects
%   <a href="matlab:help m\utilities\date2obs">m\utilities\date2obs</a>                                 - start=The_data.TimeInfo(1).date_2_observation(obj.options.estim_start_date);
%   <a href="matlab:help m\utilities\date2serial">m\utilities\date2serial</a>                              - test=(date2serial(1990):date2serial(1995))
%   <a href="matlab:help m\utilities\decipher_error">m\utilities\decipher_error</a>                           - % ====== evaluating the system ====== %
%   <a href="matlab:help m\utilities\estimation_engine">m\utilities\estimation_engine</a>                        - pt.optimset=optimset('Display','iter',...%[ off | iter | iter-detailed | notify | notify-detailed | final | final-detailed ]
%   m\utilities\find_nearest                             - (No help available)
%   <a href="matlab:help m\utilities\greek_symbols">m\utilities\greek_symbols</a>                            - %
%   <a href="matlab:help m\utilities\haver2rise">m\utilities\haver2rise</a>                               - % haver2rise.m
%   m\utilities\is_one                                   - (No help available)
%   m\utilities\is_zero                                  - (No help available)
%   m\utilities\ivech                                    - (No help available)
%   <a href="matlab:help m\utilities\list_opened_files">m\utilities\list_opened_files</a>                        - lists all the files currently opened in the editor
%   <a href="matlab:help m\utilities\locate_permutation">m\utilities\locate_permutation</a>                       - % [111 211 311 121 221 321 131 231 331 112 212 312 122 222 322 132 232
%   <a href="matlab:help m\utilities\locate_variables">m\utilities\locate_variables</a>                         - % I remove spaces in the variables just to make sure... I hope this
%   <a href="matlab:help m\utilities\mergestructures">m\utilities\mergestructures</a>                          -  structures with unique fields.
%   <a href="matlab:help m\utilities\mygrid">m\utilities\mygrid</a>                                   - v is a vector of the number of states in each dimension
%   m\utilities\mypermutation                            - (No help available)
%   m\utilities\mysetfield                               - (No help available)
%   <a href="matlab:help m\utilities\object2cell">m\utilities\object2cell</a>                              - nobj=numel(obj);
%   m\utilities\obs2date                                 - (No help available)
%   <a href="matlab:help m\utilities\online_function_evaluator">m\utilities\online_function_evaluator</a>                - evaluates string code as a function with inputs and outputs
%   <a href="matlab:help m\utilities\par_save">m\utilities\par_save</a>                                 - unction par_save(filename,variables,variables_names) %#ok<INUSL>
%   <a href="matlab:help m\utilities\parfor_save">m\utilities\parfor_save</a>                              - unction parfor_save(filename,x) %#ok<INUSD>
%   m\utilities\preserve                                 - (No help available)
%   <a href="matlab:help m\utilities\randword">m\utilities\randword</a>                                 - ymbols={'_'}; %'§','@','£','#','$','&',
%   m\utilities\resort                                   - (No help available)
%   <a href="matlab:help m\utilities\rise_saveaspdf">m\utilities\rise_saveaspdf</a>                           - unction correction=rise_saveaspdf(fig,filename)%%
%   m\utilities\save_objects                             - (No help available)
%   <a href="matlab:help m\utilities\serial2date">m\utilities\serial2date</a>                              - test=(date2serial(1990):date2serial(1995))
%   <a href="matlab:help m\utilities\simulate_shocks">m\utilities\simulate_shocks</a>                          -  shocks and keeps the effects separated. 
%   <a href="matlab:help m\utilities\time_frequency_stamp">m\utilities\time_frequency_stamp</a>                     - stamps helps put a stamp on the serial numbers such that the frequency is
%   <a href="matlab:help m\utilities\var2vma">m\utilities\var2vma</a>                                  - Computes the Vector Moving Average representation of a reduced-form VAR
%   <a href="matlab:help m\utilities\var_likelihood">m\utilities\var_likelihood</a>                           - Vi=V\eye(n);
%   m\utilities\var_ols                                  - (No help available)
%   m\utilities\vech                                     - (No help available)
%   <a href="matlab:help m\utilities\write_function_to_disk">m\utilities\write_function_to_disk</a>                   - outstring=[outstring,',res',int2str(ii)]; %#ok<AGROW>
%   <a href="matlab:help m\utilities\xl2databases">m\utilities\xl2databases</a>                             - {
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\utilities\functional_programming   %%%%%%%%%%%%%%%%%%%%
%
%   m\utilities\functional_programming\if_elseif    - (No help available)
%   m\utilities\functional_programming\if_then_else - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\utilities\kronecker_operations   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\utilities\kronecker_operations\A_kronecker_B_times_x">m\utilities\kronecker_operations\A_kronecker_B_times_x</a>        - computes the kron(A,B)*x
%   <a href="matlab:help m\utilities\kronecker_operations\A_times_k_kron_B">m\utilities\kronecker_operations\A_times_k_kron_B</a>             - computes C=A*(BxBx...xB)
%   <a href="matlab:help m\utilities\kronecker_operations\A_times_kron_B_I">m\utilities\kronecker_operations\A_times_kron_B_I</a>             - computes C=A*kron(B,Iq)
%   <a href="matlab:help m\utilities\kronecker_operations\A_times_kron_I_B">m\utilities\kronecker_operations\A_times_kron_I_B</a>             - computes C=A*kron(Iq,B)
%   <a href="matlab:help m\utilities\kronecker_operations\A_times_kron_I_B_I">m\utilities\kronecker_operations\A_times_kron_I_B_I</a>           - computes C=A*kron(kron(Iq,B),Ir)
%   <a href="matlab:help m\utilities\kronecker_operations\korder_matrix_vector">m\utilities\kronecker_operations\korder_matrix_vector</a>         - =reshape(x,[ru,cu,h]); % rename U to x in order to save on memory...
%   <a href="matlab:help m\utilities\kronecker_operations\kronecker_times_vector">m\utilities\kronecker_operations\kronecker_times_vector</a>       - this function computes kron(T1,T2)*vec(X)
%   m\utilities\kronecker_operations\kronecker_times_vector_tests - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\utilities\plottingTools   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\utilities\plottingTools\dim">m\utilities\plottingTools\dim</a>                                  - function []=shade(start,finish,colorstr);
%   <a href="matlab:help m\utilities\plottingTools\dim_relevant">m\utilities\plottingTools\dim_relevant</a>                         - shadenber.m
%   <a href="matlab:help m\utilities\plottingTools\multiWaitbar">m\utilities\plottingTools\multiWaitbar</a>                         - multiWaitbar: add, remove or update an entry on the multi waitbar
%   <a href="matlab:help m\utilities\plottingTools\nber_dates">m\utilities\plottingTools\nber_dates</a>                           - Official NBER recession dates, from
%   m\utilities\plottingTools\number_of_rows_and_columns_in_figure - (No help available)
%   <a href="matlab:help m\utilities\plottingTools\plot_decomp">m\utilities\plottingTools\plot_decomp</a>                          - ar(this_pos,'stack') %  area(this_pos)
%   <a href="matlab:help m\utilities\plottingTools\plot_packages">m\utilities\plottingTools\plot_packages</a>                        - str=['(',sprintf('%0.0f',fig),')'];
%   <a href="matlab:help m\utilities\plottingTools\plot_real_time">m\utilities\plottingTools\plot_real_time</a>                       - lot(pp.xdatenums(1:nr)',dd(:,1),'linewidth',2)%,'color',[0,0,0]
%   <a href="matlab:help m\utilities\plottingTools\plot_specs">m\utilities\plottingTools\plot_specs</a>                           - see also plot_real_time
%   <a href="matlab:help m\utilities\plottingTools\rotateXLabels">m\utilities\plottingTools\rotateXLabels</a>                        - rotateXLabels: rotate any xticklabels
%   <a href="matlab:help m\utilities\plottingTools\sup_label">m\utilities\plottingTools\sup_label</a>                            - This modifies Ben Barrowes' suplabel (see Matlab Central).
%   m\utilities\plottingTools\xrotate                              - (No help available)
%
%
