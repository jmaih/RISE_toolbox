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
