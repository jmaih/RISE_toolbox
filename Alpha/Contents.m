% RISE Toolbox
% Version alpha_20130907 (R2010b) 07-Sept-2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toolbox        RISE
%
% Version        alpha_20130907 07-Sept-2013
%
% Contents path  C:\Users\Junior\Dropbox\RISE\RISE_Toolbox
%
%
%%%%%%%%%%%%%%%%%%%%   path:    %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help rise_startup">rise_startup</a> - % this function sets up RISE (it replaces setpaths)
%   <a href="matlab:help setpaths">setpaths</a>     - if flag % then remove the paths
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
%   <a href="matlab:help classes\+distributions\laplace">classes\+distributions\laplace</a>                  - unction varargout=laplace(plb,pub,prob,c,d)% double exponential
%   <a href="matlab:help classes\+distributions\logistic">classes\+distributions\logistic</a>                 - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\lognormal">classes\+distributions\lognormal</a>                - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\normal">classes\+distributions\normal</a>                   - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\pareto">classes\+distributions\pareto</a>                   - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\truncated_normal">classes\+distributions\truncated_normal</a>         - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\uniform">classes\+distributions\uniform</a>                  - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\weibull">classes\+distributions\weibull</a>                  - The problem to solve is the following:
%   <a href="matlab:help classes\+distributions\wishart">classes\+distributions\wishart</a>                  - u=0; % <--- zeros(k,v)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+distributions\+old   %%%%%%%%%%%%%%%%%%%%
%
%   classes\+distributions\+old\encode_distribution      - (No help available)
%   <a href="matlab:help classes\+distributions\+old\inv_gammacdf">classes\+distributions\+old\inv_gammacdf</a>             - This function computes the cumulative density for the inverse gamma
%   classes\+distributions\+old\log_prior_density        - (No help available)
%   classes\+distributions\+old\log_truncated_cauchy_pdf - (No help available)
%   <a href="matlab:help classes\+distributions\+old\parameter_bounds">classes\+distributions\+old\parameter_bounds</a>         - adapted from dynare's prior_bounds
%   <a href="matlab:help classes\+distributions\+old\prior_density">classes\+distributions\+old\prior_density</a>            - % normally the encoding is not needed here
%   classes\+distributions\+old\truncated_cauchy_cdf     - (No help available)
%   classes\+distributions\+old\truncated_cauchy_inv     - (No help available)
%   classes\+distributions\+old\truncated_cauchy_pdf     - (No help available)
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
%   <a href="matlab:help classes\+msre_solvers\functional_iteration">classes\+msre_solvers\functional_iteration</a>        - based on the system X+inv(Aplus*X+A0)*Aminus=0. This
%   <a href="matlab:help classes\+msre_solvers\functional_iteration_">classes\+msre_solvers\functional_iteration_</a>       - based on the system X+inv(Aplus*X+A0)*Aminus=0. This
%   <a href="matlab:help classes\+msre_solvers\fwz_newton_system">classes\+msre_solvers\fwz_newton_system</a>           - We need to rewrite the system in the form
%   <a href="matlab:help classes\+msre_solvers\newton_kronecker">classes\+msre_solvers\newton_kronecker</a>            - alternative newton algorithm based on the expansion of the system
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
%%%%%%%%%%%%%%%%%%%%   path: classes\+ols   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+ols\ordinary_least_squares">classes\+ols\ordinary_least_squares</a> - ig2=RSS/T; % maximum likelihood estimator
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\+quasi_monte_carlo   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\+quasi_monte_carlo\halton">classes\+quasi_monte_carlo\halton</a>          - Examples:
%   classes\+quasi_monte_carlo\latin_hypercube - (No help available)
%   <a href="matlab:help classes\+quasi_monte_carlo\sobol">classes\+quasi_monte_carlo\sobol</a>           - Examples:
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@hdmr   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@hdmr\hdmr">classes\@hdmr\hdmr</a> - % objective is either : f and theta or a function that will help
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@hdmr\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@hdmr\private\orthonormal_polynomial">classes\@hdmr\private\orthonormal_polynomial</a>     - arning('off')%#ok<WNOFF> %,'MATLAB:divideByZero'
%   <a href="matlab:help classes\@hdmr\private\orthonormal_polynomial_old">classes\@hdmr\private\orthonormal_polynomial_old</a> - ==================
%   <a href="matlab:help classes\@hdmr\private\theta_to_x">classes\@hdmr\private\theta_to_x</a>                 - normalize
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise\assign_estimates">classes\@rise\assign_estimates</a>                         - this is the routine called during estimation for assigning new estimates
%   <a href="matlab:help classes\@rise\check_optimum">classes\@rise\check_optimum</a>                            - unction check_optimum(obj,varlist)%,r0,c0,N
%   <a href="matlab:help classes\@rise\counterfactual">classes\@rise\counterfactual</a>                           - shocks_db is a rise_time_series object with alternative history for the shocks
%   classes\@rise\draw_parameter                           - (No help available)
%   <a href="matlab:help classes\@rise\estimate">classes\@rise\estimate</a>                                 - varargin: data,optimset,optimizer,hessian_type,estim_start_from_mode,estim_parallel
%   <a href="matlab:help classes\@rise\evaluate">classes\@rise\evaluate</a>                                 - % re-initialize the object
%   <a href="matlab:help classes\@rise\filter">classes\@rise\filter</a>                                   - solve the object
%   <a href="matlab:help classes\@rise\forecast">classes\@rise\forecast</a>                                 - unction [cond_fkst_db,cond_fkst_mean_db,uncond_fkst_db]=forecast(obj,varargin)% historical_db
%   <a href="matlab:help classes\@rise\historical_decomposition">classes\@rise\historical_decomposition</a>                 - PURPOSE: Computes historical decompositions of a DSGE model
%   <a href="matlab:help classes\@rise\irf">classes\@rise\irf</a>                                      - 2- If the model is estimated, one may want to draw from the distribution,
%   <a href="matlab:help classes\@rise\is_stable_system">classes\@rise\is_stable_system</a>                         - this function checks that the solved system is
%   <a href="matlab:help classes\@rise\load_parameters">classes\@rise\load_parameters</a>                          - This loads the parameters and writes them to the parameter
%   <a href="matlab:help classes\@rise\log_marginal_data_density_chib_jeliazkov">classes\@rise\log_marginal_data_density_chib_jeliazkov</a> - now sample J parameter vectors from the proposal density. I don't think I
%   <a href="matlab:help classes\@rise\log_marginal_data_density_mhm">classes\@rise\log_marginal_data_density_mhm</a>            - arning off %#ok<WNOFF>
%   <a href="matlab:help classes\@rise\log_posterior_kernel">classes\@rise\log_posterior_kernel</a>                     - estim_hyperparams=obj.estim_hyperparams;
%   <a href="matlab:help classes\@rise\log_prior_density">classes\@rise\log_prior_density</a>                        - % alternatives are
%   <a href="matlab:help classes\@rise\monte_carlo_filtering">classes\@rise\monte_carlo_filtering</a>                    - pval_cutoff=0.1;
%   <a href="matlab:help classes\@rise\posterior_marginal_and_prior_densities">classes\@rise\posterior_marginal_and_prior_densities</a>   - =================
%   <a href="matlab:help classes\@rise\posterior_simulator">classes\@rise\posterior_simulator</a>                      - for ii=1:length(opt)
%   <a href="matlab:help classes\@rise\print_estimates">classes\@rise\print_estimates</a>                          - params=strcat(param_names,' : ',param_vals); This will not preserve
%   <a href="matlab:help classes\@rise\print_estimation_results">classes\@rise\print_estimation_results</a>                 - recision='%8.4f';
%   <a href="matlab:help classes\@rise\print_solution">classes\@rise\print_solution</a>                           - recision='%8.6f';
%   <a href="matlab:help classes\@rise\prior_plots">classes\@rise\prior_plots</a>                              - aveUnderName0=cell(1,nobj);%
%   <a href="matlab:help classes\@rise\rise">classes\@rise\rise</a>                                     - Description
%   <a href="matlab:help classes\@rise\set_options">classes\@rise\set_options</a>                              - function obj=set_options(obj,varargin)
%   <a href="matlab:help classes\@rise\set_parameters">classes\@rise\set_parameters</a>                           - There are 5 ways to assign parameters, depending on what exactly you
%   classes\@rise\set_properties                           - (No help available)
%   <a href="matlab:help classes\@rise\simulate">classes\@rise\simulate</a>                                 - UNTITLED6 Summary of this function goes here
%   <a href="matlab:help classes\@rise\simulation_diagnostics">classes\@rise\simulation_diagnostics</a>                   - first determine the number of chains and the number  of matrices in each
%   <a href="matlab:help classes\@rise\solve">classes\@rise\solve</a>                                    - solver options are 0 or 'msre_klein' (for constant-parameter models)
%   <a href="matlab:help classes\@rise\solve_alternatives">classes\@rise\solve_alternatives</a>                       - syntax is allobj=solve_alternatives(obj,solver,file2save2)
%   <a href="matlab:help classes\@rise\stoch_simul">classes\@rise\stoch_simul</a>                              - should also allow for passing options
%   <a href="matlab:help classes\@rise\theoretical_autocorrelations">classes\@rise\theoretical_autocorrelations</a>             - unction [A,info]=theoretical_autocorrelations(obj,ar)%,resolve_flag
%   <a href="matlab:help classes\@rise\theoretical_autocovariances">classes\@rise\theoretical_autocovariances</a>              - unction [A,retcode]=theoretical_autocovariances(obj,ar)%,resolve_flag
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise\private   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@rise\private\alpha_probability                                      - (No help available)
%   <a href="matlab:help classes\@rise\private\draw_parameters">classes\@rise\private\draw_parameters</a>                                        - this function draws a parameter and assigns it.
%   <a href="matlab:help classes\@rise\private\dsge_var_irf">classes\@rise\private\dsge_var_irf</a>                                           - compute PHIb, SIGb and ZZi
%   <a href="matlab:help classes\@rise\private\format_parameters">classes\@rise\private\format_parameters</a>                                      - re-order the names of the markov chains right when creating it...
%   <a href="matlab:help classes\@rise\private\generate_starting_point">classes\@rise\private\generate_starting_point</a>                                - % this function attempts to reduce the number of rejection of randomly
%   <a href="matlab:help classes\@rise\private\load_data">classes\@rise\private\load_data</a>                                              - unction [obj,issue,retcode]=load_data(obj,varargin)%,estimation_flag
%   <a href="matlab:help classes\@rise\private\load_functions">classes\@rise\private\load_functions</a>                                         - % % % % % % % % if ~exist([obj.options.results_folder,filesep,'macros'],'dir')
%   <a href="matlab:help classes\@rise\private\load_mode">classes\@rise\private\load_mode</a>                                              - % read the mode file
%   <a href="matlab:help classes\@rise\private\parameters_posterior_moments">classes\@rise\private\parameters_posterior_moments</a>                           - quantiles
%   classes\@rise\private\potential_scale_reduction                              - (No help available)
%   <a href="matlab:help classes\@rise\private\recast_loose_commitment_solution_into_markov_switching">classes\@rise\private\recast_loose_commitment_solution_into_markov_switching</a> - this function puts the loose commitment solution into a markov switching
%   <a href="matlab:help classes\@rise\private\save_filters - Copy">classes\@rise\private\save_filters - Copy</a>                                    - startdate=obj.options.estim_start_date; % obj.options.data.start;
%   <a href="matlab:help classes\@rise\private\save_filters">classes\@rise\private\save_filters</a>                                           - if ismember(p,[nobs,nobs+1]) % this is not very robust...
%   classes\@rise\private\store_probabilities                                    - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_anonymous   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_anonymous\rise_anonymous">classes\@rise_anonymous\rise_anonymous</a> - val_i=eval(obj(ii),varargin{:}); %#ok<*EVLC>
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_date   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_date\rise_date">classes\@rise_date\rise_date</a> - % TODO:
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_endo_priors   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_endo_priors\rise_endo_priors">classes\@rise_endo_priors\rise_endo_priors</a> - Shat % zero-frequency spectral density
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_equation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_equation\rise_equation">classes\@rise_equation\rise_equation</a> - % set property utility
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_estim_param   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_estim_param\plot">classes\@rise_estim_param\plot</a>             - % functions of the distribution
%   <a href="matlab:help classes\@rise_estim_param\rise_estim_param">classes\@rise_estim_param\rise_estim_param</a> - % %         function_log_pdf
%   classes\@rise_estim_param\set_properties   - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_param   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_param\rise_param">classes\@rise_param\rise_param</a>     - % set property utility
%   classes\@rise_param\set_properties - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_sad   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_sad\diff">classes\@rise_sad\diff</a>               - rise_sad/DIFF overloads diff with a rise_sad object argument
%   classes\@rise_sad\hessian            - (No help available)
%   classes\@rise_sad\jacobian           - (No help available)
%   <a href="matlab:help classes\@rise_sad\optimize">classes\@rise_sad\optimize</a>           - replaces
%   <a href="matlab:help classes\@rise_sad\rise_sad - Copy">classes\@rise_sad\rise_sad - Copy</a>    - % Symbolic automatic differentiation
%   <a href="matlab:help classes\@rise_sad\rise_sad">classes\@rise_sad\rise_sad</a>           - % Symbolic automatic differentiation
%   <a href="matlab:help classes\@rise_sad\rise_sad_no_delete">classes\@rise_sad\rise_sad_no_delete</a> - % Symbolic automatic differentiation
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_sad\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_sad\private\is_atom">classes\@rise_sad\private\is_atom</a>        - atom:= does not contain (+-*/^,)
%   <a href="matlab:help classes\@rise_sad\private\tryevaluate">classes\@rise_sad\private\tryevaluate</a>    - checks whether a string can be evaluated
%   <a href="matlab:help classes\@rise_sad\private\valid_varnames">classes\@rise_sad\private\valid_varnames</a> - if any(strcmp(str(1:under_score-1),vnames)) % <--- ismember(str(1:under_score-1),{'x','y','p','ss'})
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_time_series   %%%%%%%%%%%%%%%%%%%%
%
%   classes\@rise_time_series\aggregate                 - (No help available)
%   <a href="matlab:help classes\@rise_time_series\automatic_model_selection">classes\@rise_time_series\automatic_model_selection</a> - unction FinalResults=automatic_model_selection(obj,endo_name,exo_names,options)%[final_mod,x]
%   <a href="matlab:help classes\@rise_time_series\line">classes\@rise_time_series\line</a>                      - 10             'yyyy'                   2000
%   classes\@rise_time_series\mpower                    - (No help available)
%   <a href="matlab:help classes\@rise_time_series\plot">classes\@rise_time_series\plot</a>                      - 10             'yyyy'                   2000
%   <a href="matlab:help classes\@rise_time_series\plot_separate">classes\@rise_time_series\plot_separate</a>             - obj=this.subsref(S); % <--- obj=this(vnames{id}); does not work
%   <a href="matlab:help classes\@rise_time_series\rise_time_series">classes\@rise_time_series\rise_time_series</a>          - for kk=1:numel(s)
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@rise_variable   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@rise_variable\rise_variable">classes\@rise_variable\rise_variable</a> - % value moved under protection because it is necessary to pass it
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@sad   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@sad\diff">classes\@sad\diff</a>     - sad/DIFF overloads diff with a sad object argument
%   classes\@sad\hessian  - (No help available)
%   classes\@sad\jacobian - (No help available)
%   <a href="matlab:help classes\@sad\sad">classes\@sad\sad</a>      - record_book={iter,var,expression}
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@sad\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@sad\private\is_atom">classes\@sad\private\is_atom</a>        - atom:= does not contain (+-*/^,)
%   <a href="matlab:help classes\@sad\private\tryevaluate">classes\@sad\private\tryevaluate</a>    - checks whether a string can be evaluated
%   <a href="matlab:help classes\@sad\private\valid_varnames">classes\@sad\private\valid_varnames</a> - if any(strcmp(str(1:under_score-1),vnames)) % <--- ismember(str(1:under_score-1),{'x','y','p','ss'})
%
%
%%%%%%%%%%%%%%%%%%%%   path: classes\@sad_tree   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help classes\@sad_tree\Jacobian">classes\@sad_tree\Jacobian</a> - myfunc='exp(x(1)+acos(x(2)*log(x(3)))+atan(x(1)*x(2)))';
%   <a href="matlab:help classes\@sad_tree\char__">classes\@sad_tree\char__</a>   -  itself is already taken care of
%   <a href="matlab:help classes\@sad_tree\sad_tree">classes\@sad_tree\sad_tree</a> - % % % %         varagout=char(varargin)
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
%%%%%%%%%%%%%%%%%%%%   path: examples\MarkovSwitching\FarmerWaggonerZha2010   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\MarkovSwitching\FarmerWaggonerZha2010\fwz_test_suite_dynare_conf">examples\MarkovSwitching\FarmerWaggonerZha2010\fwz_test_suite_dynare_conf</a>                              - % housekeeping
%   examples\MarkovSwitching\FarmerWaggonerZha2010\generate_markov_switching_rational_expectations_problem - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\FarmerWaggonerZha2010\howto">examples\MarkovSwitching\FarmerWaggonerZha2010\howto</a>                                                   - % housekeeping
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
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_dynamic">examples\MarkovSwitching\RamirezWaggonerZha2012\frwz_dynamic</a>                                - Automagically generated on 08-Oct-2012 21:41:09
%   examples\MarkovSwitching\RamirezWaggonerZha2012\functional_iteration_convergence_conditions - (No help available)
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\howto_frwz_nk">examples\MarkovSwitching\RamirezWaggonerZha2012\howto_frwz_nk</a>                               - % Housekeeping
%   <a href="matlab:help examples\MarkovSwitching\RamirezWaggonerZha2012\rwz_steady_state">examples\MarkovSwitching\RamirezWaggonerZha2012\rwz_steady_state</a>                            - unction [ss,retcode,imposed]=rwz_steady_state(param_obj,flag)%param_struct=struct(parlist,parvals)
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
%   <a href="matlab:help examples\ModelsWithSteadyStateFile\steady_state_4_Canonical_Const">examples\ModelsWithSteadyStateFile\steady_state_4_Canonical_Const</a> - unction [ss,retcode]=steady_state_4_Canonical_Const(param_obj,flag)%param_struct=struct(parlist,parvals)
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
%%%%%%%%%%%%%%%%%%%%   path: examples\StickyInformation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\StickyInformation\howto">examples\StickyInformation\howto</a> - % housekeeping
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\StochasticReplanning_switching_Estimation   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\howto">examples\StochasticReplanning_switching_Estimation\howto</a>               - % housekeeping
%   <a href="matlab:help examples\StochasticReplanning_switching_Estimation\usmodel_steadystate">examples\StochasticReplanning_switching_Estimation\usmodel_steadystate</a> - computes the steady state for the observed variables in the smets-wouters
%
%
%%%%%%%%%%%%%%%%%%%%   path: examples\TimeSeriesObjects   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help examples\TimeSeriesObjects\howto">examples\TimeSeriesObjects\howto</a> - % create empty time series
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
%%%%%%%%%%%%%%%%%%%%   path: m\differentiation\automatic   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\differentiation\automatic\automatic">m\differentiation\automatic\automatic</a>                         - % inspired from autodiff (Matlab central)
%   <a href="matlab:help m\differentiation\automatic\derivative_indices">m\differentiation\automatic\derivative_indices</a>                - this function returns the indices for the derivatives calculated by
%   <a href="matlab:help m\differentiation\automatic\dyn_series">m\differentiation\automatic\dyn_series</a>                        - %DYN_SERIES dyn_series class constructor
%   <a href="matlab:help m\differentiation\automatic\first_order_derivatives">m\differentiation\automatic\first_order_derivatives</a>           - % activate variable ii
%   m\differentiation\automatic\more_test                         - (No help available)
%   <a href="matlab:help m\differentiation\automatic\multivariate_taylor_approximation">m\differentiation\automatic\multivariate_taylor_approximation</a> - MVT main program: returns a cell array of all multivariate Taylor Coefficients
%   <a href="matlab:help m\differentiation\automatic\semi_automatic">m\differentiation\automatic\semi_automatic</a>                    - function [x,dx] = semi_automatic_get(ad,mode)
%   <a href="matlab:help m\differentiation\automatic\test_system_fun">m\differentiation\automatic\test_system_fun</a>                   - y10=x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+c(10)*x(10);
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\differentiation\symbolic   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\differentiation\symbolic\symbolic_derivatives">m\differentiation\symbolic\symbolic_derivatives</a>               - the dictionary contains all the names of the parameters, the endogenous,
%   <a href="matlab:help m\differentiation\symbolic\symbolic_derivatives_2">m\differentiation\symbolic\symbolic_derivatives_2</a>             - 1- replace leads and lags by something else
%   <a href="matlab:help m\differentiation\symbolic\symbolic_derivator">m\differentiation\symbolic\symbolic_derivator</a>                 - the dictionary contains all the names of the parameters, the endogenous,
%   <a href="matlab:help m\differentiation\symbolic\symbolic_derivator_up_to_3rd_order">m\differentiation\symbolic\symbolic_derivator_up_to_3rd_order</a> - the dictionary contains all the names of the parameters, the endogenous,
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\filtering   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\filtering\CheckPositiveDefiniteness">m\filtering\CheckPositiveDefiniteness</a>                - all( all( M == M' ) ) & min( eig( M ) ) > 0
%   m\filtering\conditional_likelihood                   - (No help available)
%   m\filtering\initial_markov_distribution              - (No help available)
%   <a href="matlab:help m\filtering\kalman_initialization">m\filtering\kalman_initialization</a>                    - There is no documentation of this function yet.
%   <a href="matlab:help m\filtering\kalman_prediction">m\filtering\kalman_prediction</a>                        - only compute the places where there is some action
%   <a href="matlab:help m\filtering\kalman_update">m\filtering\kalman_update</a>                            - K=P(:,obs_id);
%   <a href="matlab:help m\filtering\likelihood_dsge_var">m\filtering\likelihood_dsge_var</a>                      - N.B: This function assumes the likelihood is to be maximized
%   <a href="matlab:help m\filtering\likelihood_markov_switching_dsge">m\filtering\likelihood_markov_switching_dsge</a>         - unction [LogLik,Incr,retcode,obj]=likelihood_markov_switching_dsge(params,obj)%
%   <a href="matlab:help m\filtering\likelihood_optimal_simple_rule">m\filtering\likelihood_optimal_simple_rule</a>           - % this important output is not created yet
%   <a href="matlab:help m\filtering\markov_switching_kalman_filter">m\filtering\markov_switching_kalman_filter</a>           - Detailed explanation to come here
%   <a href="matlab:help m\filtering\markov_switching_kalman_filter_real_time">m\filtering\markov_switching_kalman_filter_real_time</a> - y,... % data
%   <a href="matlab:help m\filtering\smoothing_step">m\filtering\smoothing_step</a>                           - Note that Durbin and Koopman define K=T*P*Z'*iF, while here it is defined
%   <a href="matlab:help m\filtering\symmetrize">m\filtering\symmetrize</a>                               - unction B=symmetrize(A)%,flag
%   <a href="matlab:help m\filtering\update_and_collapse">m\filtering\update_and_collapse</a>                      - Ptt(:,:,snow)=Ptt(:,:,snow)/PAItt(snow); %symmetrize()
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
%   <a href="matlab:help m\optimizers\bbo_3">m\optimizers\bbo_3</a>                        - 'mutation_type','gaussian' % 'cauchy','uniform'
%   <a href="matlab:help m\optimizers\bbo_gate">m\optimizers\bbo_gate</a>                     - Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   <a href="matlab:help m\optimizers\bbo_gate_conclude">m\optimizers\bbo_gate_conclude</a>            - exitflag=0; % disp('Too many function evaluations or iterations.')
%   <a href="matlab:help m\optimizers\bee">m\optimizers\bee</a>                          - % optimizer-specific properties
%   <a href="matlab:help m\optimizers\bee_2">m\optimizers\bee_2</a>                        - Reference: Inspired from Karaboga
%   <a href="matlab:help m\optimizers\bee_demo">m\optimizers\bee_demo</a>                     - % optimizer-specific properties
%   <a href="matlab:help m\optimizers\bee_demo_launch">m\optimizers\bee_demo_launch</a>              - =2; % number of parameters
%   <a href="matlab:help m\optimizers\bee_gate">m\optimizers\bee_gate</a>                     -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\bee_gate_conclude">m\optimizers\bee_gate_conclude</a>            - exitflag=0; % disp('Too many function evaluations or iterations.')
%   <a href="matlab:help m\optimizers\cmsa">m\optimizers\cmsa</a>                         - % algorithm specific options
%   <a href="matlab:help m\optimizers\cmsa_gate">m\optimizers\cmsa_gate</a>                    -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\cmsa_gate_conclude">m\optimizers\cmsa_gate_conclude</a>           - exitflag=0; % disp('Too many function evaluations or iterations.')
%   <a href="matlab:help m\optimizers\gampc">m\optimizers\gampc</a>                        - % algorithm specific options
%   <a href="matlab:help m\optimizers\gampc_2">m\optimizers\gampc_2</a>                      -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\hybrid_artificial_bee_colony">m\optimizers\hybrid_artificial_bee_colony</a> - 'F',0.5 % F=2*rand; %[0,2]
%   <a href="matlab:help m\optimizers\local_optimize">m\optimizers\local_optimize</a>               - trim everything
%   <a href="matlab:help m\optimizers\local_optimize_2">m\optimizers\local_optimize_2</a>             - trim everything
%   <a href="matlab:help m\optimizers\local_optimize_gate">m\optimizers\local_optimize_gate</a>          - if obj.iter>=obj.max_iter || ...
%   <a href="matlab:help m\optimizers\local_reoptimize">m\optimizers\local_reoptimize</a>             - trim everything
%   <a href="matlab:help m\optimizers\optimtestfun">m\optimizers\optimtestfun</a>                 - function u=u_func(x,a,k,m)
%   m\optimizers\rebuild_population           - (No help available)
%   <a href="matlab:help m\optimizers\studga">m\optimizers\studga</a>                       - % optimizer-specific options
%   <a href="matlab:help m\optimizers\studga_gate">m\optimizers\studga_gate</a>                  -  attempts to find the global minimum of a constrained function of
%   <a href="matlab:help m\optimizers\studga_gate_conclude">m\optimizers\studga_gate_conclude</a>         - exitflag=0; % disp('Too many function evaluations or iterations.')
%   m\optimizers\studga_gate_iris             - (No help available)
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\optimizers\private   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\optimizers\private\check_convergence">m\optimizers\private\check_convergence</a>   - if rem(obj.iter,100)==0 || ~isempty(stopflag)
%   <a href="matlab:help m\optimizers\private\clear_duplicates">m\optimizers\private\clear_duplicates</a>    - % select the chromosomes to change randomly
%   m\optimizers\private\compute_fitness     - (No help available)
%   m\optimizers\private\dispersion          - (No help available)
%   <a href="matlab:help m\optimizers\private\display_progress">m\optimizers\private\display_progress</a>    - fprintf(1,'restart # %3.0f   iter: %6.0f   fmin(global) %8.4f    fmin(iter) %8.4f    dispersion %8.4f    fcount  %8.0f   routine %s\n',...
%   m\optimizers\private\distance            - (No help available)
%   m\optimizers\private\find_farthest       - (No help available)
%   m\optimizers\private\find_nearest        - (No help available)
%   <a href="matlab:help m\optimizers\private\generate_candidates">m\optimizers\private\generate_candidates</a> -  we get 4 outputs?
%   <a href="matlab:help m\optimizers\private\manual_stopping">m\optimizers\private\manual_stopping</a>     - rawfile = char(textread(ManualStoppingFile,'%s','delimiter','\n','whitespace','','bufsize',40000));
%   m\optimizers\private\recenter            - (No help available)
%   <a href="matlab:help m\optimizers\private\uniform_sampling">m\optimizers\private\uniform_sampling</a>    - samples without repetition
%   <a href="matlab:help m\optimizers\private\weighted_sampling">m\optimizers\private\weighted_sampling</a>   - samples without repetition
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\parsers   %%%%%%%%%%%%%%%%%%%%
%
%   <a href="matlab:help m\parsers\CompileModelFile">m\parsers\CompileModelFile</a>          - by the way, I can still declare exogenous and make them observable at the
%   <a href="matlab:help m\parsers\analytical_symbolic_form">m\parsers\analytical_symbolic_form</a>  - turns x(1),y(1),ss(1),p(1) into x_1,y_1,ss_1,p_1
%   <a href="matlab:help m\parsers\cellstr2mat">m\parsers\cellstr2mat</a>               - unction out=cellstr2mat(CellMat) %,MatName
%   <a href="matlab:help m\parsers\chain_grid">m\parsers\chain_grid</a>                - v is a vector of the number of states in each markov chain
%   m\parsers\is_transition_probability - (No help available)
%   <a href="matlab:help m\parsers\symbolic2model">m\parsers\symbolic2model</a>            - shad=[shad,eqtn]; %#ok<*AGROW>
%   <a href="matlab:help m\parsers\update_markov_chains_info">m\parsers\update_markov_chains_info</a> - creates a markov chain info and updates it as new information comes in.
%
%
%%%%%%%%%%%%%%%%%%%%   path: m\solvers   %%%%%%%%%%%%%%%%%%%%
%
%   m\solvers\collapse_array                               - (No help available)
%   <a href="matlab:help m\solvers\computational_savings - Copy">m\solvers\computational_savings - Copy</a>                 - this function separates static variables from dynamic ones. It places all
%   <a href="matlab:help m\solvers\computational_savings">m\solvers\computational_savings</a>                        - this function separates static variables from dynamic ones. It places all
%   <a href="matlab:help m\solvers\dsge_lc_solve">m\solvers\dsge_lc_solve</a>                                - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\dsge_solve_aim">m\solvers\dsge_solve_aim</a>                               - ags=1; % no of lags and leads
%   <a href="matlab:help m\solvers\dsge_solve_gensys">m\solvers\dsge_solve_gensys</a>                            - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\dsge_solve_klein">m\solvers\dsge_solve_klein</a>                             - this function solves the rational expectations model
%   <a href="matlab:help m\solvers\expand_array">m\solvers\expand_array</a>                                 - case 1	% A0, Aminus or Aplus
%   <a href="matlab:help m\solvers\fix_point_iterator">m\solvers\fix_point_iterator</a>                           - this function solves for a fix point. Inputs are:
%   <a href="matlab:help m\solvers\gensys">m\solvers\gensys</a>                                       - function [G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,pi,div)
%   <a href="matlab:help m\solvers\get_default_optimization_option">m\solvers\get_default_optimization_option</a>              - reconvexify:
%   m\solvers\is_eigenvalue_solver_candidate               - (No help available)
%   <a href="matlab:help m\solvers\loose_commitment_solver">m\solvers\loose_commitment_solver</a>                      - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   <a href="matlab:help m\solvers\loose_commitment_solver_fix_point_unfinished">m\solvers\loose_commitment_solver_fix_point_unfinished</a> - Reference: Debortoli, Maih, Nunes (2010): "Loose Commitment in Medium
%   m\solvers\markov_switching_dsge_objective              - (No help available)
%   <a href="matlab:help m\solvers\markov_switching_dsge_stack">m\solvers\markov_switching_dsge_stack</a>                  - C_st,... % endo_nbr x h matrix of constant
%   m\solvers\msre_aim                                     - (No help available)
%   <a href="matlab:help m\solvers\msre_gensys">m\solvers\msre_gensys</a>                                  - if ~isempty(B)
%   <a href="matlab:help m\solvers\msre_klein">m\solvers\msre_klein</a>                                   - if ~isempty(B)
%   <a href="matlab:help m\solvers\msre_matrix_times_vector">m\solvers\msre_matrix_times_vector</a>                     - the old version is faster and solves but solves the problem
%   <a href="matlab:help m\solvers\msre_solve">m\solvers\msre_solve</a>                                   - This procedure assumes the steady state has been solved and that apart
%   <a href="matlab:help m\solvers\msre_solve_accelerated">m\solvers\msre_solve_accelerated</a>                       - This procedure assumes the steady state has been solved and that apart
%   <a href="matlab:help m\solvers\qzdiv">m\solvers\qzdiv</a>                                        - function [A,B,Q,Z] = qzdiv(stake,A,B,Q,Z)
%   <a href="matlab:help m\solvers\qzswitch">m\solvers\qzswitch</a>                                     - function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%   <a href="matlab:help m\solvers\schur_solver">m\solvers\schur_solver</a>                                 - T,S,Q,Z] = qz(F,G,'real');%complex
%   <a href="matlab:help m\solvers\solve_steady_state">m\solvers\solve_steady_state</a>                           - % compute the constant
%   <a href="matlab:help m\solvers\transpose_free_quasi_minimum_residual">m\solvers\transpose_free_quasi_minimum_residual</a>        - A,... % coefficient matrix
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
%   <a href="matlab:help m\utilities\A_kronecker_B_times_x">m\utilities\A_kronecker_B_times_x</a>                - computes the kron(A,B)*x
%   m\utilities\CheckArgument                        - (No help available)
%   <a href="matlab:help m\utilities\cell2object">m\utilities\cell2object</a>                          - ypes=C(3,:); %#ok<NASGU>
%   <a href="matlab:help m\utilities\decipher_error">m\utilities\decipher_error</a>                       - % ====== evaluating the system ====== %
%   <a href="matlab:help m\utilities\dim">m\utilities\dim</a>                                  - function []=shade(start,finish,colorstr);
%   <a href="matlab:help m\utilities\dim_relevant">m\utilities\dim_relevant</a>                         - shadenber.m
%   <a href="matlab:help m\utilities\estimation_engine">m\utilities\estimation_engine</a>                    - pt.optimset=optimset('Display','iter',...%[ off | iter | iter-detailed | notify | notify-detailed | final | final-detailed ]
%   m\utilities\find_nearest                         - (No help available)
%   <a href="matlab:help m\utilities\greek_symbols">m\utilities\greek_symbols</a>                        - %
%   <a href="matlab:help m\utilities\haver2rise">m\utilities\haver2rise</a>                           - % haver2rise.m
%   m\utilities\ivech                                - (No help available)
%   <a href="matlab:help m\utilities\kronecker_times_vector">m\utilities\kronecker_times_vector</a>               - this function computes kron(T1,T2)*vec(X)
%   m\utilities\kronecker_times_vector_tests         - (No help available)
%   <a href="matlab:help m\utilities\locate_variables">m\utilities\locate_variables</a>                     - % I remove spaces in the variables just to make sure... I hope this
%   <a href="matlab:help m\utilities\mergestructures">m\utilities\mergestructures</a>                      -  structures with unique fields.
%   <a href="matlab:help m\utilities\mygrid">m\utilities\mygrid</a>                               - v is a vector of the number of states in each dimension
%   m\utilities\mysetfield                           - (No help available)
%   <a href="matlab:help m\utilities\nber_dates">m\utilities\nber_dates</a>                           - Official NBER recession dates, from
%   m\utilities\number_of_rows_and_columns_in_figure - (No help available)
%   <a href="matlab:help m\utilities\object2cell">m\utilities\object2cell</a>                          - nobj=numel(obj);
%   <a href="matlab:help m\utilities\par_save">m\utilities\par_save</a>                             - unction par_save(filename,variables,variables_names) %#ok<INUSL>
%   <a href="matlab:help m\utilities\parfor_save">m\utilities\parfor_save</a>                          - unction parfor_save(filename,x) %#ok<INUSD>
%   m\utilities\preserve                             - (No help available)
%   <a href="matlab:help m\utilities\randword">m\utilities\randword</a>                             - ymbols={'_'}; %'§','@','£','#','$','&',
%   <a href="matlab:help m\utilities\recursive_moments">m\utilities\recursive_moments</a>                    - {
%   m\utilities\resort                               - (No help available)
%   m\utilities\save_objects                         - (No help available)
%   <a href="matlab:help m\utilities\sup_label">m\utilities\sup_label</a>                            - This modifies Ben Barrowes' suplabel (see Matlab Central).
%   <a href="matlab:help m\utilities\test_is_stable_system_old">m\utilities\test_is_stable_system_old</a>            - this function checks that the solved system is
%   <a href="matlab:help m\utilities\var_likelihood">m\utilities\var_likelihood</a>                       - Vi=V\eye(n);
%   m\utilities\var_ols                              - (No help available)
%   m\utilities\vech                                 - (No help available)
%   <a href="matlab:help m\utilities\write_function_to_disk">m\utilities\write_function_to_disk</a>               - outstring=[outstring,',res',int2str(ii)]; %#ok<AGROW>
%   <a href="matlab:help m\utilities\xl2databases">m\utilities\xl2databases</a>                         - {
%
%
