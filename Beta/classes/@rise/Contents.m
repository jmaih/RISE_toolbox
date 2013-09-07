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
