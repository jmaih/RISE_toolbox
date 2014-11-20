% RISE Toolbox
% Version 1.0.1 (R2014b) 20-Nov-2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Toolbox        RISE
%
% Version        1.0.1 20-Nov-2014
%
% Contents path  C:\Users\jma\Documents\GitHub\RISE_toolbox\help
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox   -------------------%
%
%   <a href="matlab:help rise_exit">rise_exit</a>                                              -  exit from RISE
%   <a href="matlab:help rise_startup">rise_startup</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\log_marginal_data_density   -------------------%
%
%   <a href="matlab:help classes\log_marginal_data_density\chib_jeliazkov">classes\log_marginal_data_density\chib_jeliazkov</a>                                                   -  computation of the MDD using chib and jeliazkov
%   <a href="matlab:help classes\log_marginal_data_density\modified_harmonic_mean">classes\log_marginal_data_density\modified_harmonic_mean</a>                                           -  computation of the MDD using modified harmonic mean
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\+conditional_forecasting   -------------------%
%
%   <a href="matlab:help classes\models\+conditional_forecasting\parse_information">classes\models\+conditional_forecasting\parse_information</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\+dsge_tools   -------------------%
%
%   <a href="matlab:help classes\models\+dsge_tools\ergodic_parameters">classes\models\+dsge_tools\ergodic_parameters</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\+dsge_tools\+utils   -------------------%
%
%   <a href="matlab:help classes\models\+dsge_tools\+utils\msre_initial_guess">classes\models\+dsge_tools\+utils\msre_initial_guess</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\+generic_tools   -------------------%
%
%   <a href="matlab:help classes\models\+generic_tools\choose_state">classes\models\+generic_tools\choose_state</a>                                                             - H1 line
%   <a href="matlab:help classes\models\+generic_tools\parameter_position_and_regimes">classes\models\+generic_tools\parameter_position_and_regimes</a>                                           - H1 line
%   <a href="matlab:help classes\models\+generic_tools\set_exogenous_data">classes\models\+generic_tools\set_exogenous_data</a>                                                       - H1 line
%   <a href="matlab:help classes\models\+generic_tools\set_simulation_regimes">classes\models\+generic_tools\set_simulation_regimes</a>                                                   - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\+vartools   -------------------%
%
%   <a href="matlab:help classes\models\+vartools\build_parameter_vector">classes\models\+vartools\build_parameter_vector</a>                                                      - H1 line
%   <a href="matlab:help classes\models\+vartools\companion_form">classes\models\+vartools\companion_form</a>                                                              - H1 line
%   <a href="matlab:help classes\models\+vartools\covariance_decomposition">classes\models\+vartools\covariance_decomposition</a>                                                    - H1 line
%   <a href="matlab:help classes\models\+vartools\create_baseline_parameters">classes\models\+vartools\create_baseline_parameters</a>                                                  - H1 line
%   <a href="matlab:help classes\models\+vartools\decompose_impact">classes\models\+vartools\decompose_impact</a>                                                            - H1 line
%   <a href="matlab:help classes\models\+vartools\find_long_run">classes\models\+vartools\find_long_run</a>                                                               - H1 line
%   <a href="matlab:help classes\models\+vartools\format_markov_chains">classes\models\+vartools\format_markov_chains</a>                                                        - H1 line
%   <a href="matlab:help classes\models\+vartools\ivech">classes\models\+vartools\ivech</a>                                                                       - H1 line
%   <a href="matlab:help classes\models\+vartools\moving_average_representation">classes\models\+vartools\moving_average_representation</a>                                               - H1 line
%   <a href="matlab:help classes\models\+vartools\ols">classes\models\+vartools\ols</a>                                                                         - H1 line
%   <a href="matlab:help classes\models\+vartools\parameters_to_matrices">classes\models\+vartools\parameters_to_matrices</a>                                                      - H1 line
%   <a href="matlab:help classes\models\+vartools\preliminaries">classes\models\+vartools\preliminaries</a>                                                               - H1 line
%   <a href="matlab:help classes\models\+vartools\reset_markov_chains">classes\models\+vartools\reset_markov_chains</a>                                                         - H1 line
%   <a href="matlab:help classes\models\+vartools\resolve">classes\models\+vartools\resolve</a>                                                                     - H1 line
%   <a href="matlab:help classes\models\+vartools\select_parameter_type">classes\models\+vartools\select_parameter_type</a>                                                       - H1 line
%   <a href="matlab:help classes\models\+vartools\set_structural_matrices_structure">classes\models\+vartools\set_structural_matrices_structure</a>                                           - H1 line
%   <a href="matlab:help classes\models\+vartools\set_y_and_x">classes\models\+vartools\set_y_and_x</a>                                                                 - H1 line
%   <a href="matlab:help classes\models\+vartools\var_likelihood">classes\models\+vartools\var_likelihood</a>                                                              - H1 line
%   <a href="matlab:help classes\models\+vartools\vec">classes\models\+vartools\vec</a>                                                                         - H1 line
%   <a href="matlab:help classes\models\+vartools\vech">classes\models\+vartools\vech</a>                                                                        - H1 line
%   <a href="matlab:help classes\models\+vartools\vma">classes\models\+vartools\vma</a>                                                                         - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\dsge   -------------------%
%
%   <a href="matlab:help classes\models\dsge\check_derivatives">classes\models\dsge\check_derivatives</a>                                                          -  - compares the derivatives and the solutions from various differentiation techniques
%   <a href="matlab:help classes\models\dsge\compute_steady_state">classes\models\dsge\compute_steady_state</a>                                                       - H1 line
%   classes\models\dsge\conclude_estimation                                                        - (No help available)
%   <a href="matlab:help classes\models\dsge\create_estimation_blocks">classes\models\dsge\create_estimation_blocks</a>                                                   - H1 line
%   <a href="matlab:help classes\models\dsge\create_state_list">classes\models\dsge\create_state_list</a>                                                          -  creates the list of the state variables in the solution
%   classes\models\dsge\do_not_anticipate_future_shocks                                            - (No help available)
%   classes\models\dsge\dsge                                                                       - (No help available)
%   <a href="matlab:help classes\models\dsge\filter">classes\models\dsge\filter</a>                                                                     - H1 line
%   <a href="matlab:help classes\models\dsge\forecast_real_time">classes\models\dsge\forecast_real_time</a>                                                         -  - forecast from each point in time
%   <a href="matlab:help classes\models\dsge\irf">classes\models\dsge\irf</a>                                                                        - H1 line
%   <a href="matlab:help classes\models\dsge\is_stable_system">classes\models\dsge\is_stable_system</a>                                                           - H1 line
%   classes\models\dsge\latex_model_file                                                           - (No help available)
%   classes\models\dsge\load_data                                                                  - (No help available)
%   classes\models\dsge\load_solution                                                              - (No help available)
%   <a href="matlab:help classes\models\dsge\monte_carlo_filtering">classes\models\dsge\monte_carlo_filtering</a>                                                      - H1 line
%   classes\models\dsge\prepare_transition_routine                                                 - (No help available)
%   <a href="matlab:help classes\models\dsge\print_solution">classes\models\dsge\print_solution</a>                                                             -  - print the solution of a model or vector of models
%   classes\models\dsge\re_order_output_rows                                                       - (No help available)
%   <a href="matlab:help classes\models\dsge\resid">classes\models\dsge\resid</a>                                                                      - H1 line
%   <a href="matlab:help classes\models\dsge\set">classes\models\dsge\set</a>                                                                        -  - sets options for dsge|rise models
%   <a href="matlab:help classes\models\dsge\set_solution_to_companion">classes\models\dsge\set_solution_to_companion</a>                                                  - H1 line
%   classes\models\dsge\set_z_eplus_horizon                                                        - (No help available)
%   classes\models\dsge\simulate                                                                   - (No help available)
%   <a href="matlab:help classes\models\dsge\simulate_nonlinear">classes\models\dsge\simulate_nonlinear</a>                                                         - H1 line
%   <a href="matlab:help classes\models\dsge\solve">classes\models\dsge\solve</a>                                                                      - H1 line
%   <a href="matlab:help classes\models\dsge\solve_alternatives">classes\models\dsge\solve_alternatives</a>                                                         - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\favar   -------------------%
%
%   <a href="matlab:help classes\models\favar\conclude_estimation">classes\models\favar\conclude_estimation</a>                                            - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\rfvar   -------------------%
%
%   <a href="matlab:help classes\models\rfvar\check_identification">classes\models\rfvar\check_identification</a>                                                                      - H1 line
%   classes\models\rfvar\find_posterior_mode                                                                       - (No help available)
%   classes\models\rfvar\initialize_posterior_simulation                                                           - (No help available)
%   classes\models\rfvar\rfvar                                                                                     - (No help available)
%   classes\models\rfvar\set_structural_shocks                                                                     - (No help available)
%   <a href="matlab:help classes\models\rfvar\solve">classes\models\rfvar\solve</a>                                                                                     - H1 line
%   classes\models\rfvar\sort_Q                                                                                    - (No help available)
%   <a href="matlab:help classes\models\rfvar\structural_form">classes\models\rfvar\structural_form</a>                                                                           -  finds A structural form given the imposed restrictions
%   classes\models\rfvar\translate_restrictions                                                                    - (No help available)
%   classes\models\rfvar\update_estimated_parameter_names                                                          - (No help available)
%   classes\models\rfvar\update_posterior_simulation_initial_conditions                                            - (No help available)
%   classes\models\rfvar\var_rotation                                                                              - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\rise   -------------------%
%
%   classes\models\rise\rise                                            - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\rise_generic   -------------------%
%
%   classes\models\rise_generic\add_markov_chains_and_parameters                                                          - (No help available)
%   classes\models\rise_generic\assign_estimates                                                                          - (No help available)
%   <a href="matlab:help classes\models\rise_generic\check_optimum">classes\models\rise_generic\check_optimum</a>                                                                             - H1 line
%   <a href="matlab:help classes\models\rise_generic\decompose_parameter_name">classes\models\rise_generic\decompose_parameter_name</a>                                                                  - H1 line
%   classes\models\rise_generic\do_names                                                                                  - (No help available)
%   <a href="matlab:help classes\models\rise_generic\draw_parameter">classes\models\rise_generic\draw_parameter</a>                                                                            - H1 line
%   <a href="matlab:help classes\models\rise_generic\estimate">classes\models\rise_generic\estimate</a>                                                                                  -  - estimates the parameters of a RISE model
%   classes\models\rise_generic\estimation_wrapper                                                                        - (No help available)
%   classes\models\rise_generic\find_posterior_mode                                                                       - (No help available)
%   <a href="matlab:help classes\models\rise_generic\forecast">classes\models\rise_generic\forecast</a>                                                                                  -  - computes forecasts for rise|dsge|svar|rfvar models
%   <a href="matlab:help classes\models\rise_generic\get">classes\models\rise_generic\get</a>                                                                                       - H1 line
%   classes\models\rise_generic\get_estimated_parameter_names                                                             - (No help available)
%   <a href="matlab:help classes\models\rise_generic\historical_decomposition">classes\models\rise_generic\historical_decomposition</a>                                                                  -  Computes historical decompositions of a DSGE model
%   classes\models\rise_generic\initialize_posterior_simulation                                                           - (No help available)
%   <a href="matlab:help classes\models\rise_generic\irf">classes\models\rise_generic\irf</a>                                                                                       -  - computes impulse responses for a RISE model
%   <a href="matlab:help classes\models\rise_generic\isnan">classes\models\rise_generic\isnan</a>                                                                                     - H1 line
%   classes\models\rise_generic\load_data                                                                                 - (No help available)
%   classes\models\rise_generic\load_mode                                                                                 - (No help available)
%   <a href="matlab:help classes\models\rise_generic\load_parameters">classes\models\rise_generic\load_parameters</a>                                                                           - H1 line
%   <a href="matlab:help classes\models\rise_generic\log_marginal_data_density">classes\models\rise_generic\log_marginal_data_density</a>                                                                 - H1 line
%   <a href="matlab:help classes\models\rise_generic\log_posterior_kernel">classes\models\rise_generic\log_posterior_kernel</a>                                                                      - H1 line
%   <a href="matlab:help classes\models\rise_generic\log_prior_density">classes\models\rise_generic\log_prior_density</a>                                                                         - H1 line
%   <a href="matlab:help classes\models\rise_generic\parameters_links">classes\models\rise_generic\parameters_links</a>                                                                          - H1 line
%   <a href="matlab:help classes\models\rise_generic\posterior_marginal_and_prior_densities">classes\models\rise_generic\posterior_marginal_and_prior_densities</a>                                                    - H1 line
%   <a href="matlab:help classes\models\rise_generic\posterior_simulator">classes\models\rise_generic\posterior_simulator</a>                                                                       - H1 line
%   classes\models\rise_generic\prepare_transition_routine                                                                - (No help available)
%   <a href="matlab:help classes\models\rise_generic\print_estimation_results">classes\models\rise_generic\print_estimation_results</a>                                                                  - H1 line
%   <a href="matlab:help classes\models\rise_generic\prior_plots">classes\models\rise_generic\prior_plots</a>                                                                               - H1 line
%   classes\models\rise_generic\re_order_output_rows                                                                      - (No help available)
%   <a href="matlab:help classes\models\rise_generic\refresh">classes\models\rise_generic\refresh</a>                                                                                   -  - refresh the options of an old object with a newer version of
%   <a href="matlab:help classes\models\rise_generic\report">classes\models\rise_generic\report</a>                                                                                    -  assigns the elements of interest to a rise_report.report object
%   classes\models\rise_generic\rise_generic                                                                              - (No help available)
%   <a href="matlab:help classes\models\rise_generic\set">classes\models\rise_generic\set</a>                                                                                       -  - sets options for RISE models
%   classes\models\rise_generic\set_simulation_initial_conditions                                                         - (No help available)
%   <a href="matlab:help classes\models\rise_generic\setup_calibration">classes\models\rise_generic\setup_calibration</a>                                                                         - H1 line
%   <a href="matlab:help classes\models\rise_generic\setup_general_restrictions">classes\models\rise_generic\setup_general_restrictions</a>                                                                - H1 line
%   <a href="matlab:help classes\models\rise_generic\setup_linear_restrictions">classes\models\rise_generic\setup_linear_restrictions</a>                                                                 - H1 line
%   <a href="matlab:help classes\models\rise_generic\setup_priors">classes\models\rise_generic\setup_priors</a>                                                                              - H1 line
%   <a href="matlab:help classes\models\rise_generic\simulate">classes\models\rise_generic\simulate</a>                                                                                  -  - simulates a RISE model
%   <a href="matlab:help classes\models\rise_generic\simulation_diagnostics">classes\models\rise_generic\simulation_diagnostics</a>                                                                    - H1 line
%   <a href="matlab:help classes\models\rise_generic\stoch_simul">classes\models\rise_generic\stoch_simul</a>                                                                               - H1 line
%   <a href="matlab:help classes\models\rise_generic\theoretical_autocorrelations">classes\models\rise_generic\theoretical_autocorrelations</a>                                                              - H1 line
%   <a href="matlab:help classes\models\rise_generic\theoretical_autocovariances">classes\models\rise_generic\theoretical_autocovariances</a>                                                               - H1 line
%   classes\models\rise_generic\update_estimated_parameter_names                                                          - (No help available)
%   classes\models\rise_generic\update_posterior_simulation_initial_conditions                                            - (No help available)
%   <a href="matlab:help classes\models\rise_generic\variance_decomposition">classes\models\rise_generic\variance_decomposition</a>                                                                    - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\stochvol   -------------------%
%
%   classes\models\stochvol\form_parameter_matrices                                              - (No help available)
%   classes\models\stochvol\format_blocks                                                        - (No help available)
%   classes\models\stochvol\redo_declarations                                                    - (No help available)
%   <a href="matlab:help classes\models\stochvol\setup_linear_restrictions">classes\models\stochvol\setup_linear_restrictions</a>                                            - H1 line
%   classes\models\stochvol\simulation_engine                                                    - (No help available)
%   classes\models\stochvol\stochvol                                                             - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\models\svar   -------------------%
%
%   classes\models\svar\conclude_estimation                                                         - (No help available)
%   <a href="matlab:help classes\models\svar\create_estimated_parameters_list">classes\models\svar\create_estimated_parameters_list</a>                                            - H1 line
%   classes\models\svar\decompose_parameter                                                         - (No help available)
%   classes\models\svar\load_solution                                                               - (No help available)
%   <a href="matlab:help classes\models\svar\msvar_priors">classes\models\svar\msvar_priors</a>                                                                - H1 line
%   classes\models\svar\reformat_restriction                                                        - (No help available)
%   <a href="matlab:help classes\models\svar\set_solution_to_companion">classes\models\svar\set_solution_to_companion</a>                                                   - H1 line
%   classes\models\svar\simulation_engine                                                           - (No help available)
%   <a href="matlab:help classes\models\svar\solve">classes\models\svar\solve</a>                                                                       - H1 line
%   classes\models\svar\svar                                                                        - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\time_series   -------------------%
%
%   <a href="matlab:help classes\time_series\char2serial">classes\time_series\char2serial</a>                                                    - H1 line
%   <a href="matlab:help classes\time_series\date2obs">classes\time_series\date2obs</a>                                                       - H1 line
%   <a href="matlab:help classes\time_series\date2serial">classes\time_series\date2serial</a>                                                    - H1 line
%   <a href="matlab:help classes\time_series\date2year_period">classes\time_series\date2year_period</a>                                               - H1 line
%   <a href="matlab:help classes\time_series\frequency2char">classes\time_series\frequency2char</a>                                                 - H1 line
%   <a href="matlab:help classes\time_series\frequency2num">classes\time_series\frequency2num</a>                                                  - H1 line
%   <a href="matlab:help classes\time_series\frequency_map">classes\time_series\frequency_map</a>                                                  - H1 line
%   <a href="matlab:help classes\time_series\is_serial">classes\time_series\is_serial</a>                                                      - H1 line
%   <a href="matlab:help classes\time_series\obs2date">classes\time_series\obs2date</a>                                                       - H1 line
%   <a href="matlab:help classes\time_series\parse_plot_args">classes\time_series\parse_plot_args</a>                                                - H1 line
%   <a href="matlab:help classes\time_series\period2period">classes\time_series\period2period</a>                                                  - H1 line
%   <a href="matlab:help classes\time_series\serial2date">classes\time_series\serial2date</a>                                                    - H1 line
%   <a href="matlab:help classes\time_series\serial2frequency">classes\time_series\serial2frequency</a>                                               - H1 line
%   <a href="matlab:help classes\time_series\time_frequency_stamp">classes\time_series\time_frequency_stamp</a>                                           - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\time_series\rseries   -------------------%
%
%   classes\time_series\rseries\rseries                                            - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\time_series\ts   -------------------%
%
%   <a href="matlab:help classes\time_series\ts\acos">classes\time_series\ts\acos</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\acosh">classes\time_series\ts\acosh</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\acot">classes\time_series\ts\acot</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\acoth">classes\time_series\ts\acoth</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\aggregate">classes\time_series\ts\aggregate</a>                                                            - H1 line
%   <a href="matlab:help classes\time_series\ts\allmean">classes\time_series\ts\allmean</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\and">classes\time_series\ts\and</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\apply">classes\time_series\ts\apply</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\asin">classes\time_series\ts\asin</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\asinh">classes\time_series\ts\asinh</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\atan">classes\time_series\ts\atan</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\atanh">classes\time_series\ts\atanh</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\automatic_model_selection">classes\time_series\ts\automatic_model_selection</a>                                            - H1 line
%   <a href="matlab:help classes\time_series\ts\bar">classes\time_series\ts\bar</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\barh">classes\time_series\ts\barh</a>                                                                 - H1 line
%   classes\time_series\ts\binary_operation                                                     - (No help available)
%   <a href="matlab:help classes\time_series\ts\boxplot">classes\time_series\ts\boxplot</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\bsxfun">classes\time_series\ts\bsxfun</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\cat">classes\time_series\ts\cat</a>                                                                  -  concatenates time series along the specified dimension
%   <a href="matlab:help classes\time_series\ts\collect">classes\time_series\ts\collect</a>                                                              -  - brings together several time series object into a one time series
%   classes\time_series\ts\comparison                                                           - (No help available)
%   <a href="matlab:help classes\time_series\ts\corr">classes\time_series\ts\corr</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\corrcoef">classes\time_series\ts\corrcoef</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\cos">classes\time_series\ts\cos</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\cosh">classes\time_series\ts\cosh</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\cot">classes\time_series\ts\cot</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\coth">classes\time_series\ts\coth</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\cov">classes\time_series\ts\cov</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\ctranspose">classes\time_series\ts\ctranspose</a>                                                           - H1 line
%   <a href="matlab:help classes\time_series\ts\cumprod">classes\time_series\ts\cumprod</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\cumsum">classes\time_series\ts\cumsum</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\decompose_series">classes\time_series\ts\decompose_series</a>                                                     - H1 line
%   <a href="matlab:help classes\time_series\ts\describe">classes\time_series\ts\describe</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\display">classes\time_series\ts\display</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\double">classes\time_series\ts\double</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\drop">classes\time_series\ts\drop</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\dummy">classes\time_series\ts\dummy</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\eq">classes\time_series\ts\eq</a>                                                                   - H1 line
%   <a href="matlab:help classes\time_series\ts\exp">classes\time_series\ts\exp</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\expanding">classes\time_series\ts\expanding</a>                                                            - H1 line
%   <a href="matlab:help classes\time_series\ts\fanchart">classes\time_series\ts\fanchart</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\ge">classes\time_series\ts\ge</a>                                                                   - H1 line
%   <a href="matlab:help classes\time_series\ts\get">classes\time_series\ts\get</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\gt">classes\time_series\ts\gt</a>                                                                   - H1 line
%   <a href="matlab:help classes\time_series\ts\head">classes\time_series\ts\head</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\hist">classes\time_series\ts\hist</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\horzcat">classes\time_series\ts\horzcat</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\hpfilter">classes\time_series\ts\hpfilter</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\index">classes\time_series\ts\index</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\interpolate">classes\time_series\ts\interpolate</a>                                                          - H1 line
%   <a href="matlab:help classes\time_series\ts\intersect">classes\time_series\ts\intersect</a>                                                            - H1 line
%   <a href="matlab:help classes\time_series\ts\isfinite">classes\time_series\ts\isfinite</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\isinf">classes\time_series\ts\isinf</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\isnan">classes\time_series\ts\isnan</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\jbtest">classes\time_series\ts\jbtest</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\kurtosis">classes\time_series\ts\kurtosis</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\le">classes\time_series\ts\le</a>                                                                   - H1 line
%   <a href="matlab:help classes\time_series\ts\log">classes\time_series\ts\log</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\lt">classes\time_series\ts\lt</a>                                                                   - H1 line
%   classes\time_series\ts\main_frame                                                           - (No help available)
%   <a href="matlab:help classes\time_series\ts\max">classes\time_series\ts\max</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\mean">classes\time_series\ts\mean</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\median">classes\time_series\ts\median</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\min">classes\time_series\ts\min</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\minus">classes\time_series\ts\minus</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\mode">classes\time_series\ts\mode</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\mpower">classes\time_series\ts\mpower</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\mrdivide">classes\time_series\ts\mrdivide</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\mtimes">classes\time_series\ts\mtimes</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\nan">classes\time_series\ts\nan</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\ne">classes\time_series\ts\ne</a>                                                                   - H1 line
%   <a href="matlab:help classes\time_series\ts\numel">classes\time_series\ts\numel</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\ones">classes\time_series\ts\ones</a>                                                                 -  overloads ones for ts objects
%   <a href="matlab:help classes\time_series\ts\pages2struct">classes\time_series\ts\pages2struct</a>                                                         - H1 line
%   <a href="matlab:help classes\time_series\ts\plot">classes\time_series\ts\plot</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\plotyy">classes\time_series\ts\plotyy</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\plus">classes\time_series\ts\plus</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\power">classes\time_series\ts\power</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\prctile">classes\time_series\ts\prctile</a>                                                              -  Percentiles of a time series (ts)
%   classes\time_series\ts\process_subs                                                         - (No help available)
%   <a href="matlab:help classes\time_series\ts\quantile">classes\time_series\ts\quantile</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\rand">classes\time_series\ts\rand</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\randn">classes\time_series\ts\randn</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\range">classes\time_series\ts\range</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\rdivide">classes\time_series\ts\rdivide</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\regress">classes\time_series\ts\regress</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\reset_start_date">classes\time_series\ts\reset_start_date</a>                                                     - H1 line
%   <a href="matlab:help classes\time_series\ts\rolling">classes\time_series\ts\rolling</a>                                                              - H1 line
%   classes\time_series\ts\set_locations                                                        - (No help available)
%   <a href="matlab:help classes\time_series\ts\sin">classes\time_series\ts\sin</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\sinh">classes\time_series\ts\sinh</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\skewness">classes\time_series\ts\skewness</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\sort">classes\time_series\ts\sort</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\spectrum">classes\time_series\ts\spectrum</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\std">classes\time_series\ts\std</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\step_dummy">classes\time_series\ts\step_dummy</a>                                                           - H1 line
%   <a href="matlab:help classes\time_series\ts\subsasgn">classes\time_series\ts\subsasgn</a>                                                             - H1 line
%   <a href="matlab:help classes\time_series\ts\subsref">classes\time_series\ts\subsref</a>                                                              - H1 line
%   <a href="matlab:help classes\time_series\ts\sum">classes\time_series\ts\sum</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\tail">classes\time_series\ts\tail</a>                                                                 - H1 line
%   <a href="matlab:help classes\time_series\ts\times">classes\time_series\ts\times</a>                                                                - H1 line
%   <a href="matlab:help classes\time_series\ts\transform">classes\time_series\ts\transform</a>                                                            - H1 line
%   <a href="matlab:help classes\time_series\ts\transpose">classes\time_series\ts\transpose</a>                                                            - H1 line
%   classes\time_series\ts\ts                                                                   - (No help available)
%   classes\time_series\ts\ts_roll_or_expand                                                    - (No help available)
%   <a href="matlab:help classes\time_series\ts\uminus">classes\time_series\ts\uminus</a>                                                               - H1 line
%   classes\time_series\ts\unary_operation                                                      - (No help available)
%   <a href="matlab:help classes\time_series\ts\values">classes\time_series\ts\values</a>                                                               - H1 line
%   <a href="matlab:help classes\time_series\ts\var">classes\time_series\ts\var</a>                                                                  - H1 line
%   <a href="matlab:help classes\time_series\ts\zeros">classes\time_series\ts\zeros</a>                                                                - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\utils\+quasi_monte_carlo   -------------------%
%
%   <a href="matlab:help classes\utils\+quasi_monte_carlo\halton">classes\utils\+quasi_monte_carlo\halton</a>                                                    - H1 line
%   <a href="matlab:help classes\utils\+quasi_monte_carlo\latin_hypercube">classes\utils\+quasi_monte_carlo\latin_hypercube</a>                                           - H1 line
%   <a href="matlab:help classes\utils\+quasi_monte_carlo\sobol">classes\utils\+quasi_monte_carlo\sobol</a>                                                     - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\utils\+rise_report   -------------------%
%
%   <a href="matlab:help classes\utils\+rise_report\chapter">classes\utils\+rise_report\chapter</a>                                                   -  report chapter object
%   <a href="matlab:help classes\utils\+rise_report\cleardoublepage">classes\utils\+rise_report\cleardoublepage</a>                                           -  report clear double page object
%   <a href="matlab:help classes\utils\+rise_report\clearpage">classes\utils\+rise_report\clearpage</a>                                                 -  report clear page object
%   <a href="matlab:help classes\utils\+rise_report\enumerate">classes\utils\+rise_report\enumerate</a>                                                 -  report enumerate object
%   <a href="matlab:help classes\utils\+rise_report\feed_properties">classes\utils\+rise_report\feed_properties</a>                                           - H1 line
%   <a href="matlab:help classes\utils\+rise_report\figure">classes\utils\+rise_report\figure</a>                                                    -  report figure object
%   <a href="matlab:help classes\utils\+rise_report\footnote">classes\utils\+rise_report\footnote</a>                                                  -  report footnote object
%   <a href="matlab:help classes\utils\+rise_report\generic_report">classes\utils\+rise_report\generic_report</a>                                            -  generic class for report objects
%   <a href="matlab:help classes\utils\+rise_report\include">classes\utils\+rise_report\include</a>                                                   -  report include object
%   <a href="matlab:help classes\utils\+rise_report\itemize">classes\utils\+rise_report\itemize</a>                                                   -  report itemize object
%   <a href="matlab:help classes\utils\+rise_report\newpage">classes\utils\+rise_report\newpage</a>                                                   -  report new page object
%   <a href="matlab:help classes\utils\+rise_report\page_styles">classes\utils\+rise_report\page_styles</a>                                               - H1 line
%   <a href="matlab:help classes\utils\+rise_report\pagebreak">classes\utils\+rise_report\pagebreak</a>                                                 -  report page break object
%   <a href="matlab:help classes\utils\+rise_report\paragraph">classes\utils\+rise_report\paragraph</a>                                                 -  report paragraph object
%   <a href="matlab:help classes\utils\+rise_report\quotation">classes\utils\+rise_report\quotation</a>                                                 -  report quotation object
%   <a href="matlab:help classes\utils\+rise_report\report">classes\utils\+rise_report\report</a>                                                    -  Reporting system
%   <a href="matlab:help classes\utils\+rise_report\section">classes\utils\+rise_report\section</a>                                                   -  reporting section
%   <a href="matlab:help classes\utils\+rise_report\subparagraph">classes\utils\+rise_report\subparagraph</a>                                              -  report subparagraph object
%   <a href="matlab:help classes\utils\+rise_report\subsection">classes\utils\+rise_report\subsection</a>                                                - subection report section object
%   <a href="matlab:help classes\utils\+rise_report\subsubsection">classes\utils\+rise_report\subsubsection</a>                                             -  report subsubsection object
%   <a href="matlab:help classes\utils\+rise_report\table">classes\utils\+rise_report\table</a>                                                     -  report table object
%   <a href="matlab:help classes\utils\+rise_report\text">classes\utils\+rise_report\text</a>                                                      -  report text object
%   <a href="matlab:help classes\utils\+rise_report\title_item">classes\utils\+rise_report\title_item</a>                                                - H1 line
%   <a href="matlab:help classes\utils\+rise_report\titlepage">classes\utils\+rise_report\titlepage</a>                                                 -  report title page object
%   <a href="matlab:help classes\utils\+rise_report\verbatim">classes\utils\+rise_report\verbatim</a>                                                  -  report verbatim object
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\utils\coef   -------------------%
%
%   classes\utils\coef\coef                                            - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\utils\hdmr   -------------------%
%
%   classes\utils\hdmr\hdmr                                            - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\classes\utils\mcf   -------------------%
%
%   classes\utils\mcf\mcf                                            - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+ar1_approximation   -------------------%
%
%   <a href="matlab:help m\+ar1_approximation\rouwenhorst">m\+ar1_approximation\rouwenhorst</a>                                           -  approximates and AR(1) process with a Markov chain
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+distributions   -------------------%
%
%   <a href="matlab:help m\+distributions\beta">m\+distributions\beta</a>                                                               - H1 line
%   <a href="matlab:help m\+distributions\cauchy">m\+distributions\cauchy</a>                                                             - H1 line
%   <a href="matlab:help m\+distributions\distribution_tests">m\+distributions\distribution_tests</a>                                                 -  of distributions
%   <a href="matlab:help m\+distributions\empirical_cdf">m\+distributions\empirical_cdf</a>                                                      - H1 line
%   <a href="matlab:help m\+distributions\empirical_moments">m\+distributions\empirical_moments</a>                                                  - H1 line
%   <a href="matlab:help m\+distributions\find_bounds">m\+distributions\find_bounds</a>                                                        - H1 line
%   <a href="matlab:help m\+distributions\find_hyperparameters">m\+distributions\find_hyperparameters</a>                                               - H1 line
%   <a href="matlab:help m\+distributions\gamma">m\+distributions\gamma</a>                                                              - H1 line
%   <a href="matlab:help m\+distributions\hyperparameter_residuals">m\+distributions\hyperparameter_residuals</a>                                           - H1 line
%   <a href="matlab:help m\+distributions\inv_gamma">m\+distributions\inv_gamma</a>                                                          - H1 line
%   <a href="matlab:help m\+distributions\inv_wishart">m\+distributions\inv_wishart</a>                                                        - H1 line
%   <a href="matlab:help m\+distributions\kernel_density">m\+distributions\kernel_density</a>                                                     - H1 line
%   <a href="matlab:help m\+distributions\laplace">m\+distributions\laplace</a>                                                            - H1 line
%   <a href="matlab:help m\+distributions\left_triang">m\+distributions\left_triang</a>                                                        - H1 line
%   <a href="matlab:help m\+distributions\logistic">m\+distributions\logistic</a>                                                           - H1 line
%   <a href="matlab:help m\+distributions\lognormal">m\+distributions\lognormal</a>                                                          - H1 line
%   <a href="matlab:help m\+distributions\normal">m\+distributions\normal</a>                                                             - H1 line
%   <a href="matlab:help m\+distributions\pareto">m\+distributions\pareto</a>                                                             - H1 line
%   <a href="matlab:help m\+distributions\right_triang">m\+distributions\right_triang</a>                                                       - H1 line
%   <a href="matlab:help m\+distributions\truncated_normal">m\+distributions\truncated_normal</a>                                                   - H1 line
%   <a href="matlab:help m\+distributions\uniform">m\+distributions\uniform</a>                                                            - H1 line
%   <a href="matlab:help m\+distributions\weibull">m\+distributions\weibull</a>                                                            - H1 line
%   <a href="matlab:help m\+distributions\wishart">m\+distributions\wishart</a>                                                            - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+msre_solvers   -------------------%
%
%   <a href="matlab:help m\+msre_solvers\functional_iteration_h">m\+msre_solvers\functional_iteration_h</a>                                                - H1 line
%   <a href="matlab:help m\+msre_solvers\functional_iteration_h_full">m\+msre_solvers\functional_iteration_h_full</a>                                           - H1 line
%   <a href="matlab:help m\+msre_solvers\fwz_newton_system">m\+msre_solvers\fwz_newton_system</a>                                                     - H1 line
%   <a href="matlab:help m\+msre_solvers\newton_iteration_h">m\+msre_solvers\newton_iteration_h</a>                                                    - H1 line
%   <a href="matlab:help m\+msre_solvers\newton_iteration_h_full">m\+msre_solvers\newton_iteration_h_full</a>                                               - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+optimization   -------------------%
%
%   <a href="matlab:help m\+optimization\estimation_engine">m\+optimization\estimation_engine</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+truncated_multivariate_normal   -------------------%
%
%   <a href="matlab:help m\+truncated_multivariate_normal\geweke_hajivassiliou_keane">m\+truncated_multivariate_normal\geweke_hajivassiliou_keane</a>                                           - H1 line
%   <a href="matlab:help m\+truncated_multivariate_normal\gibbs">m\+truncated_multivariate_normal\gibbs</a>                                                                - H1 line
%   <a href="matlab:help m\+truncated_multivariate_normal\quick_and_dirty">m\+truncated_multivariate_normal\quick_and_dirty</a>                                                      - H1 line
%   <a href="matlab:help m\+truncated_multivariate_normal\univariate_draw">m\+truncated_multivariate_normal\univariate_draw</a>                                                      - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+unit_tests   -------------------%
%
%   <a href="matlab:help m\+unit_tests\hessian_computation_test">m\+unit_tests\hessian_computation_test</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils   -------------------%
%
%   <a href="matlab:help m\+utils\saveaspdf">m\+utils\saveaspdf</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+code   -------------------%
%
%   <a href="matlab:help m\+utils\+code\code2file">m\+utils\+code\code2file</a>                                                                               - H1 line
%   <a href="matlab:help m\+utils\+code\code2func">m\+utils\+code\code2func</a>                                                                               - H1 line
%   <a href="matlab:help m\+utils\+code\code2vector">m\+utils\+code\code2vector</a>                                                                             -  - transforms a set of function handles into a single function
%   <a href="matlab:help m\+utils\+code\evaluate_automatic_derivatives">m\+utils\+code\evaluate_automatic_derivatives</a>                                                          - H1 line
%   <a href="matlab:help m\+utils\+code\evaluate_functions">m\+utils\+code\evaluate_functions</a>                                                                      -  - evaluates functions in various formats
%   <a href="matlab:help m\+utils\+code\evaluate_jacobian_numerically">m\+utils\+code\evaluate_jacobian_numerically</a>                                                           -  - numerical evaluation of the jacobian of the objective function
%   <a href="matlab:help m\+utils\+code\evaluate_policy_objective_hessian_numerically">m\+utils\+code\evaluate_policy_objective_hessian_numerically</a>                                           -  - numerical evaluation of the hessian of the policy objective
%   <a href="matlab:help m\+utils\+code\evaluate_transition_matrices">m\+utils\+code\evaluate_transition_matrices</a>                                                            - H1 line
%   <a href="matlab:help m\+utils\+code\func2fhandle">m\+utils\+code\func2fhandle</a>                                                                            -  - creates a handle to an m-file not on the matlab search path
%   <a href="matlab:help m\+utils\+code\validate_transition_matrix">m\+utils\+code\validate_transition_matrix</a>                                                              - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+cov   -------------------%
%
%   <a href="matlab:help m\+utils\+cov\nearest">m\+utils\+cov\nearest</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+cov\symmetrize">m\+utils\+cov\symmetrize</a>                                           -  - makes a square matrix symmetric
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+error   -------------------%
%
%   <a href="matlab:help m\+utils\+error\decipher">m\+utils\+error\decipher</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+error\valid">m\+utils\+error\valid</a>                                              - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+estim   -------------------%
%
%   <a href="matlab:help m\+utils\+estim\generate_starting_point">m\+utils\+estim\generate_starting_point</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+filtering   -------------------%
%
%   <a href="matlab:help m\+utils\+filtering\check_steady_state_kalman">m\+utils\+filtering\check_steady_state_kalman</a>                                           -  - checks whether the covariance matrix has converged
%   <a href="matlab:help m\+utils\+filtering\prediction_step">m\+utils\+filtering\prediction_step</a>                                                     -  - prediction step for Kalman filter
%   <a href="matlab:help m\+utils\+filtering\smoothing_step">m\+utils\+filtering\smoothing_step</a>                                                      -  - smoothing step
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+forecast   -------------------%
%
%   <a href="matlab:help m\+utils\+forecast\aggregate_initial_conditions">m\+utils\+forecast\aggregate_initial_conditions</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+forecast\check_shock_id">m\+utils\+forecast\check_shock_id</a>                                                         - H1 line
%   <a href="matlab:help m\+utils\+forecast\create_shocks">m\+utils\+forecast\create_shocks</a>                                                          - H1 line
%   <a href="matlab:help m\+utils\+forecast\irf">m\+utils\+forecast\irf</a>                                                                    - H1 line
%   <a href="matlab:help m\+utils\+forecast\load_start_values">m\+utils\+forecast\load_start_values</a>                                                      -  - load the start values for forecasting
%   <a href="matlab:help m\+utils\+forecast\multi_step">m\+utils\+forecast\multi_step</a>                                                             - H1 line
%   <a href="matlab:help m\+utils\+forecast\nullify_deterministic_shocks">m\+utils\+forecast\nullify_deterministic_shocks</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+forecast\one_step">m\+utils\+forecast\one_step</a>                                                               - H1 line
%   <a href="matlab:help m\+utils\+forecast\replace_impulse">m\+utils\+forecast\replace_impulse</a>                                                        - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+forecast\+conditional   -------------------%
%
%   <a href="matlab:help m\+utils\+forecast\+conditional\build_shock_restrictions">m\+utils\+forecast\+conditional\build_shock_restrictions</a>                                                                      -  - build restriction conditions on the shocks
%   m\+utils\+forecast\+conditional\compute_forecasts                                                                             - (No help available)
%   <a href="matlab:help m\+utils\+forecast\+conditional\conditional_distribution_of_endogenous_and_shocks">m\+utils\+forecast\+conditional\conditional_distribution_of_endogenous_and_shocks</a>                                             - =============================
%   <a href="matlab:help m\+utils\+forecast\+conditional\forecast_engine">m\+utils\+forecast\+conditional\forecast_engine</a>                                                                               -  - computes conditional forecasts for linear models
%   m\+utils\+forecast\+conditional\initialize_array                                                                              - (No help available)
%   <a href="matlab:help m\+utils\+forecast\+conditional\null_and_column_spaces">m\+utils\+forecast\+conditional\null_and_column_spaces</a>                                                                        - H1 line
%   <a href="matlab:help m\+utils\+forecast\+conditional\number_of_conditioning_shocks_periods">m\+utils\+forecast\+conditional\number_of_conditioning_shocks_periods</a>                                                         -  - length of shocks over the forecast horizon
%   <a href="matlab:help m\+utils\+forecast\+conditional\projection_sub_engine">m\+utils\+forecast\+conditional\projection_sub_engine</a>                                                                         - H1 line
%   <a href="matlab:help m\+utils\+forecast\+conditional\remove_holes">m\+utils\+forecast\+conditional\remove_holes</a>                                                                                  - H1 line
%   <a href="matlab:help m\+utils\+forecast\+conditional\state_matrices">m\+utils\+forecast\+conditional\state_matrices</a>                                                                                - H1 line
%   <a href="matlab:help m\+utils\+forecast\+conditional\truncated_mv_normal_rnd">m\+utils\+forecast\+conditional\truncated_mv_normal_rnd</a>                                                                       - x=truncated_mv_normal_rnd(mu,SIG,lb,ub)
%   m\+utils\+forecast\+conditional\unconditional_distribution_of_endogenous_and_shocks                                           - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+forecast\oldtools   -------------------%
%
%   m\+utils\+forecast\oldtools\forecast_engine                                                        - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_engine                                                     - (No help available)
%   m\+utils\+forecast\oldtools\three_pass_regression_filter                                           - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+forecast\oldtools\forecasting_tools   -------------------%
%
%   m\+utils\+forecast\oldtools\forecasting_tools\ComputeForecasts                                                                         - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\ConditionalDistributionOfEndogenousAndShocks                                             - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\ConditionalProjectionSubEngine                                                           - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\ConditionalStateMatrices                                                                 - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\InitializeArray                                                                          - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\NullAndColumnSpaces                                                                      - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\RemoveHoles                                                                              - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\TruncatedMultivariateNormalRnd                                                           - (No help available)
%   m\+utils\+forecast\oldtools\forecasting_tools\UnconditionalDistributionOfEndogenousAndShocks                                           - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+functional_programming   -------------------%
%
%   <a href="matlab:help m\+utils\+functional_programming\if_elseif">m\+utils\+functional_programming\if_elseif</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+functional_programming\if_then_else">m\+utils\+functional_programming\if_then_else</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+gridfuncs   -------------------%
%
%   <a href="matlab:help m\+utils\+gridfuncs\build_grid">m\+utils\+gridfuncs\build_grid</a>                                                   - H1 line
%   <a href="matlab:help m\+utils\+gridfuncs\chain_grid">m\+utils\+gridfuncs\chain_grid</a>                                                   - H1 line
%   <a href="matlab:help m\+utils\+gridfuncs\commutation">m\+utils\+gridfuncs\commutation</a>                                                  - H1 line
%   <a href="matlab:help m\+utils\+gridfuncs\locate_permutation">m\+utils\+gridfuncs\locate_permutation</a>                                           -  locates a permutation deriving from differentiation
%   <a href="matlab:help m\+utils\+gridfuncs\mygrid">m\+utils\+gridfuncs\mygrid</a>                                                       -  creates a grid of points
%   <a href="matlab:help m\+utils\+gridfuncs\mypermutation">m\+utils\+gridfuncs\mypermutation</a>                                                - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+hessian   -------------------%
%
%   <a href="matlab:help m\+utils\+hessian\conditioner">m\+utils\+hessian\conditioner</a>                                                  - H1 line
%   <a href="matlab:help m\+utils\+hessian\finite_differences">m\+utils\+hessian\finite_differences</a>                                           -  - computes the hessian by finite differences
%   <a href="matlab:help m\+utils\+hessian\outer_product">m\+utils\+hessian\outer_product</a>                                                - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+kronecker   -------------------%
%
%   <a href="matlab:help m\+utils\+kronecker\A_kron_B_times_x">m\+utils\+kronecker\A_kron_B_times_x</a>                                                       - H1 line
%   <a href="matlab:help m\+utils\+kronecker\A_times_B_kron_C">m\+utils\+kronecker\A_times_B_kron_C</a>                                                       - H1 line
%   <a href="matlab:help m\+utils\+kronecker\A_times_k_kron_B">m\+utils\+kronecker\A_times_k_kron_B</a>                                                       - H1 line
%   <a href="matlab:help m\+utils\+kronecker\A_times_kron_B_I">m\+utils\+kronecker\A_times_kron_B_I</a>                                                       - H1 line
%   <a href="matlab:help m\+utils\+kronecker\A_times_kron_I_B">m\+utils\+kronecker\A_times_kron_I_B</a>                                                       - H1 line
%   <a href="matlab:help m\+utils\+kronecker\A_times_kron_I_B_I">m\+utils\+kronecker\A_times_kron_I_B_I</a>                                                     - H1 line
%   <a href="matlab:help m\+utils\+kronecker\fernandez_plateau_stewart">m\+utils\+kronecker\fernandez_plateau_stewart</a>                                              -  computes kron(Q1,Q2,...,Qk)*x or x*kron(Q1,...,Qk)
%   <a href="matlab:help m\+utils\+kronecker\korder_matrix_vector">m\+utils\+kronecker\korder_matrix_vector</a>                                                   - H1 line
%   <a href="matlab:help m\+utils\+kronecker\kron_Q1_Qk_times_X">m\+utils\+kronecker\kron_Q1_Qk_times_X</a>                                                     - X_times_kron_Q1_Qk Efficient Kronecker Multiplication of kron(Q1,Q2,...,Qk)*X
%   <a href="matlab:help m\+utils\+kronecker\kron_times_vector">m\+utils\+kronecker\kron_times_vector</a>                                                      - H1 line
%   m\+utils\+kronecker\kronecker_times_vector_tests                                           - (No help available)
%   <a href="matlab:help m\+utils\+kronecker\perfect_shuffle">m\+utils\+kronecker\perfect_shuffle</a>                                                        -  produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%   <a href="matlab:help m\+utils\+kronecker\shrink_expand">m\+utils\+kronecker\shrink_expand</a>                                                          -  computes shrinking and expansion objects for the manipulation of symmetric tensors
%   <a href="matlab:help m\+utils\+kronecker\shuffle">m\+utils\+kronecker\shuffle</a>                                                                - perfect_shuffle produces the a perfect shuffle matrix that turns kron(B,C) into kron(C,B)
%   <a href="matlab:help m\+utils\+kronecker\sum">m\+utils\+kronecker\sum</a>                                                                    -  sums products of all possible permutations of kronecker products
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+latex   -------------------%
%
%   <a href="matlab:help m\+utils\+latex\pdflatex">m\+utils\+latex\pdflatex</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+mcmc   -------------------%
%
%   <a href="matlab:help m\+utils\+mcmc\alpha_probability">m\+utils\+mcmc\alpha_probability</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+mcmc\constant_bvar_sampler">m\+utils\+mcmc\constant_bvar_sampler</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+mcmc\mh_sampler">m\+utils\+mcmc\mh_sampler</a>                                                      - H1 line
%   <a href="matlab:help m\+utils\+mcmc\mh_sampling_engine">m\+utils\+mcmc\mh_sampling_engine</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+mcmc\parameters_moments">m\+utils\+mcmc\parameters_moments</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+mcmc\random_walk_mcmc">m\+utils\+mcmc\random_walk_mcmc</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+mcmc\retune_coefficient">m\+utils\+mcmc\retune_coefficient</a>                                              - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+miscellaneous   -------------------%
%
%   <a href="matlab:help m\+utils\+miscellaneous\CheckArgument">m\+utils\+miscellaneous\CheckArgument</a>                                                                                     - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\E_kron_X_Y">m\+utils\+miscellaneous\E_kron_X_Y</a>                                                                                        - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\cell2object">m\+utils\+miscellaneous\cell2object</a>                                                                                       - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\estimated_time_of_arrival">m\+utils\+miscellaneous\estimated_time_of_arrival</a>                                                                         - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\find_farthest">m\+utils\+miscellaneous\find_farthest</a>                                                                                     - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\find_nearest">m\+utils\+miscellaneous\find_nearest</a>                                                                                      - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\generate_markov_switching_rational_expectations_problem">m\+utils\+miscellaneous\generate_markov_switching_rational_expectations_problem</a>                                           - H1 line
%   m\+utils\+miscellaneous\greek_symbols                                                                                     - (No help available)
%   <a href="matlab:help m\+utils\+miscellaneous\integrate_regimes">m\+utils\+miscellaneous\integrate_regimes</a>                                                                                 - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\list_opened_files">m\+utils\+miscellaneous\list_opened_files</a>                                                                                 - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\mergestructures">m\+utils\+miscellaneous\mergestructures</a>                                                                                   - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\object2cell">m\+utils\+miscellaneous\object2cell</a>                                                                                       - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\parse_arguments">m\+utils\+miscellaneous\parse_arguments</a>                                                                                   - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\preserve">m\+utils\+miscellaneous\preserve</a>                                                                                          - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\randword">m\+utils\+miscellaneous\randword</a>                                                                                          - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\setfield">m\+utils\+miscellaneous\setfield</a>                                                                                          - H1 line
%   <a href="matlab:help m\+utils\+miscellaneous\sum_products">m\+utils\+miscellaneous\sum_products</a>                                                                                      - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+moments   -------------------%
%
%   <a href="matlab:help m\+utils\+moments\recursive">m\+utils\+moments\recursive</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+numdiff   -------------------%
%
%   <a href="matlab:help m\+utils\+numdiff\hessian">m\+utils\+numdiff\hessian</a>                                            -  - computes the hessian of a scalar function
%   <a href="matlab:help m\+utils\+numdiff\jacobian">m\+utils\+numdiff\jacobian</a>                                           -  - computes the jacobian of a function or a vector of functions
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+optim   -------------------%
%
%   <a href="matlab:help m\+utils\+optim\check_convergence">m\+utils\+optim\check_convergence</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+optim\clear_duplicates">m\+utils\+optim\clear_duplicates</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+optim\compare_individuals">m\+utils\+optim\compare_individuals</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+optim\compute_fitness">m\+utils\+optim\compute_fitness</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+optim\dispersion">m\+utils\+optim\dispersion</a>                                                    - H1 line
%   <a href="matlab:help m\+utils\+optim\display_progress">m\+utils\+optim\display_progress</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+optim\distance">m\+utils\+optim\distance</a>                                                      - H1 line
%   <a href="matlab:help m\+utils\+optim\dynamic_penalty">m\+utils\+optim\dynamic_penalty</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+optim\generate_candidates">m\+utils\+optim\generate_candidates</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+optim\manual_stopping">m\+utils\+optim\manual_stopping</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+optim\recenter">m\+utils\+optim\recenter</a>                                                      - H1 line
%   <a href="matlab:help m\+utils\+optim\resort">m\+utils\+optim\resort</a>                                                        - H1 line
%   <a href="matlab:help m\+utils\+optim\selection_process">m\+utils\+optim\selection_process</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+optim\sort_population">m\+utils\+optim\sort_population</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+optim\uniform_sampling">m\+utils\+optim\uniform_sampling</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+optim\universal_options">m\+utils\+optim\universal_options</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+optim\weighted_sampling">m\+utils\+optim\weighted_sampling</a>                                             - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+parallel   -------------------%
%
%   <a href="matlab:help m\+utils\+parallel\par_save">m\+utils\+parallel\par_save</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+parallel\parfor_save">m\+utils\+parallel\parfor_save</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+plot   -------------------%
%
%   <a href="matlab:help m\+utils\+plot\multiple">m\+utils\+plot\multiple</a>                                                          - H1 line
%   <a href="matlab:help m\+utils\+plot\myplot">m\+utils\+plot\myplot</a>                                                            - H1 line
%   <a href="matlab:help m\+utils\+plot\one_figure_rows_columns">m\+utils\+plot\one_figure_rows_columns</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+plot\prior_posterior">m\+utils\+plot\prior_posterior</a>                                                   - H1 line
%   <a href="matlab:help m\+utils\+plot\saveaspdf">m\+utils\+plot\saveaspdf</a>                                                         - H1 line
%   <a href="matlab:help m\+utils\+plot\waitbar">m\+utils\+plot\waitbar</a>                                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+smooth_approximation   -------------------%
%
%   <a href="matlab:help m\+utils\+smooth_approximation\sabs">m\+utils\+smooth_approximation\sabs</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+smooth_approximation\smax">m\+utils\+smooth_approximation\smax</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+smooth_approximation\smin">m\+utils\+smooth_approximation\smin</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+smooth_transition   -------------------%
%
%   <a href="matlab:help m\+utils\+smooth_transition\exponential">m\+utils\+smooth_transition\exponential</a>                                                     - H1 line
%   <a href="matlab:help m\+utils\+smooth_transition\logistic">m\+utils\+smooth_transition\logistic</a>                                                        - H1 line
%   <a href="matlab:help m\+utils\+smooth_transition\nth_order_logistic">m\+utils\+smooth_transition\nth_order_logistic</a>                                              - H1 line
%   <a href="matlab:help m\+utils\+smooth_transition\second_order_logistic">m\+utils\+smooth_transition\second_order_logistic</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+solve   -------------------%
%
%   <a href="matlab:help m\+utils\+solve\partition_variables">m\+utils\+solve\partition_variables</a>                                                   - H1 line
%   <a href="matlab:help m\+utils\+solve\pull_first_order_partitions">m\+utils\+solve\pull_first_order_partitions</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+solve\solution_topology">m\+utils\+solve\solution_topology</a>                                                     - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+stat   -------------------%
%
%   <a href="matlab:help m\+utils\+stat\corr">m\+utils\+stat\corr</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+stat\corrcoef">m\+utils\+stat\corrcoef</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+stat\cov">m\+utils\+stat\cov</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+stat\kurtosis">m\+utils\+stat\kurtosis</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+stat\max">m\+utils\+stat\max</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+stat\mean">m\+utils\+stat\mean</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+stat\median">m\+utils\+stat\median</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+stat\min">m\+utils\+stat\min</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+stat\mode">m\+utils\+stat\mode</a>                                               - H1 line
%   <a href="matlab:help m\+utils\+stat\nanmax">m\+utils\+stat\nanmax</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+stat\nanmean">m\+utils\+stat\nanmean</a>                                            - H1 line
%   <a href="matlab:help m\+utils\+stat\nanmin">m\+utils\+stat\nanmin</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+stat\nansum">m\+utils\+stat\nansum</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+stat\nanvar">m\+utils\+stat\nanvar</a>                                             - H1 line
%   <a href="matlab:help m\+utils\+stat\quantile">m\+utils\+stat\quantile</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+stat\skewness">m\+utils\+stat\skewness</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+stat\std">m\+utils\+stat\std</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+stat\sum">m\+utils\+stat\sum</a>                                                - H1 line
%   <a href="matlab:help m\+utils\+stat\var">m\+utils\+stat\var</a>                                                - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+struct   -------------------%
%
%   <a href="matlab:help m\+utils\+struct\pad">m\+utils\+struct\pad</a>                                           -  concatenates fields structures horizontally
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+table   -------------------%
%
%   <a href="matlab:help m\+utils\+table\concatenate">m\+utils\+table\concatenate</a>                                           - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\+utils\+time_series   -------------------%
%
%   <a href="matlab:help m\+utils\+time_series\concatenate_series_from_different_models">m\+utils\+time_series\concatenate_series_from_different_models</a>                                           - H1 line
%   <a href="matlab:help m\+utils\+time_series\data_request">m\+utils\+time_series\data_request</a>                                                                       -  - selects the data requested for estimation or for forecasting
%   <a href="matlab:help m\+utils\+time_series\haver2rise">m\+utils\+time_series\haver2rise</a>                                                                         - H1 line
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\differentiation\aplanar   -------------------%
%
%   m\differentiation\aplanar\aplanar                                            - (No help available)
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\differentiation\splanar   -------------------%
%
%   <a href="matlab:help m\differentiation\splanar\char">m\differentiation\splanar\char</a>                                                     -  - transforms splanar derivatives to strings
%   <a href="matlab:help m\differentiation\splanar\differentiate">m\differentiation\splanar\differentiate</a>                                            -  - differentiates vectors of splanar objects
%   <a href="matlab:help m\differentiation\splanar\get">m\differentiation\splanar\get</a>                                                      - H1 line
%   <a href="matlab:help m\differentiation\splanar\kron">m\differentiation\splanar\kron</a>                                                     - H1 line
%   <a href="matlab:help m\differentiation\splanar\load_varlist">m\differentiation\splanar\load_varlist</a>                                             - H1 line
%   <a href="matlab:help m\differentiation\splanar\print">m\differentiation\splanar\print</a>                                                    -  - transforms the output of splanar.differentiate into char or functions
%   <a href="matlab:help m\differentiation\splanar\set">m\differentiation\splanar\set</a>                                                      - H1 line
%   m\differentiation\splanar\splanar                                                  - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\filtering   -------------------%
%
%   <a href="matlab:help m\filtering\conditional_likelihood">m\filtering\conditional_likelihood</a>                                                     - H1 line
%   <a href="matlab:help m\filtering\initial_markov_distribution">m\filtering\initial_markov_distribution</a>                                                - H1 line
%   <a href="matlab:help m\filtering\kalman_initialization">m\filtering\kalman_initialization</a>                                                      - H1 line
%   <a href="matlab:help m\filtering\likelihood_dsge_var">m\filtering\likelihood_dsge_var</a>                                                        - H1 line
%   <a href="matlab:help m\filtering\likelihood_markov_switching_dsge">m\filtering\likelihood_markov_switching_dsge</a>                                           - H1 line
%   <a href="matlab:help m\filtering\likelihood_optimal_simple_rule">m\filtering\likelihood_optimal_simple_rule</a>                                             - H1 line
%   <a href="matlab:help m\filtering\msre_kalman">m\filtering\msre_kalman</a>                                                                - H1 line
%   <a href="matlab:help m\filtering\msre_kalman_cell">m\filtering\msre_kalman_cell</a>                                                           - H1 line
%   <a href="matlab:help m\filtering\msre_kalman_cell_real_time">m\filtering\msre_kalman_cell_real_time</a>                                                 - H1 line
%   <a href="matlab:help m\filtering\msre_linear_filter">m\filtering\msre_linear_filter</a>                                                         - H1 line
%   <a href="matlab:help m\filtering\simulation_smoother">m\filtering\simulation_smoother</a>                                                        - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\optimizers   -------------------%
%
%   <a href="matlab:help m\optimizers\bee_gate">m\optimizers\bee_gate</a>                                                         -  gateway to bee
%   <a href="matlab:help m\optimizers\blockwise_optimization">m\optimizers\blockwise_optimization</a>                                           -  optimization by blocks of parameters rather than the whole vector
%
%
%-------------------   class: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\optimizers\bee   -------------------%
%
%   m\optimizers\bee\bee                                            - (No help available)
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\plotting_tools   -------------------%
%
%   <a href="matlab:help m\plotting_tools\dim">m\plotting_tools\dim</a>                                                                            - H1 line
%   <a href="matlab:help m\plotting_tools\dim_relevant">m\plotting_tools\dim_relevant</a>                                                                   - H1 line
%   <a href="matlab:help m\plotting_tools\nber_dates">m\plotting_tools\nber_dates</a>                                                                     - H1 line
%   <a href="matlab:help m\plotting_tools\number_of_rows_and_columns_in_figure">m\plotting_tools\number_of_rows_and_columns_in_figure</a>                                           - H1 line
%   <a href="matlab:help m\plotting_tools\plot_decomp">m\plotting_tools\plot_decomp</a>                                                                    - H1 line
%   <a href="matlab:help m\plotting_tools\plot_fanchart">m\plotting_tools\plot_fanchart</a>                                                                  - H1 line
%   <a href="matlab:help m\plotting_tools\plot_real_time">m\plotting_tools\plot_real_time</a>                                                                 - H1 line
%   <a href="matlab:help m\plotting_tools\plot_specs">m\plotting_tools\plot_specs</a>                                                                     - H1 line
%   <a href="matlab:help m\plotting_tools\rotateXLabels">m\plotting_tools\rotateXLabels</a>                                                                  - H1 line
%   <a href="matlab:help m\plotting_tools\sup_label">m\plotting_tools\sup_label</a>                                                                      - H1 line
%   <a href="matlab:help m\plotting_tools\xrotate">m\plotting_tools\xrotate</a>                                                                        - H1 line
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\shortcuts   -------------------%
%
%   <a href="matlab:help m\shortcuts\decipher">m\shortcuts\decipher</a>                                                        - H1 line
%   <a href="matlab:help m\shortcuts\exponential">m\shortcuts\exponential</a>                                                     - H1 line
%   <a href="matlab:help m\shortcuts\if_elseif">m\shortcuts\if_elseif</a>                                                       - low-level function
%   <a href="matlab:help m\shortcuts\if_then_else">m\shortcuts\if_then_else</a>                                                    - low-level function
%   <a href="matlab:help m\shortcuts\ivech">m\shortcuts\ivech</a>                                                           - low-level function
%   <a href="matlab:help m\shortcuts\locate_variables">m\shortcuts\locate_variables</a>                                                -  locate variables in an array
%   <a href="matlab:help m\shortcuts\logistic">m\shortcuts\logistic</a>                                                        -  1st-order logistic function
%   <a href="matlab:help m\shortcuts\logistic2">m\shortcuts\logistic2</a>                                                       -  2nd-order logistic function
%   <a href="matlab:help m\shortcuts\logisticn">m\shortcuts\logisticn</a>                                                       -  nth-order logistic function
%   <a href="matlab:help m\shortcuts\newreport">m\shortcuts\newreport</a>                                                       -  initialize a new report object
%   <a href="matlab:help m\shortcuts\nth_order_logistic">m\shortcuts\nth_order_logistic</a>                                              - H1 line
%   <a href="matlab:help m\shortcuts\sabs">m\shortcuts\sabs</a>                                                            -  smooth approximation of abs
%   <a href="matlab:help m\shortcuts\second_order_logistic">m\shortcuts\second_order_logistic</a>                                           - H1 line
%   <a href="matlab:help m\shortcuts\smax">m\shortcuts\smax</a>                                                            -  smooth approximation of max(x,0)
%   <a href="matlab:help m\shortcuts\smin">m\shortcuts\smin</a>                                                            -  smooth approximation of min(x,0)
%   <a href="matlab:help m\shortcuts\vec">m\shortcuts\vec</a>                                                             - H1 line
%   <a href="matlab:help m\shortcuts\vech">m\shortcuts\vech</a>                                                            - low-level function
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\solvers   -------------------%
%
%   <a href="matlab:help m\solvers\compute_steady_state_transition_matrix">m\solvers\compute_steady_state_transition_matrix</a>                                           -  computes the transition matrix a
%   <a href="matlab:help m\solvers\fix_point_iterator">m\solvers\fix_point_iterator</a>                                                               -  solves the fix point of a function
%   <a href="matlab:help m\solvers\solve_steady_state">m\solvers\solve_steady_state</a>                                                               -  solves the steady state
%   <a href="matlab:help m\solvers\transpose_free_quasi_minimum_residual">m\solvers\transpose_free_quasi_minimum_residual</a>                                            -  attempts to solve Ax=b
%
%
%-------------------   path: C:\Users\jma\Documents\GitHub\RISE_toolbox\m\solvers\X_equal_A_X_B_plus_C_solvers   -------------------%
%
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\diffuse_lyapunov_equation">m\solvers\X_equal_A_X_B_plus_C_solvers\diffuse_lyapunov_equation</a>                                           -  attempts to solve V=T*V*T'+RQR under
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\doubling_solve">m\solvers\X_equal_A_X_B_plus_C_solvers\doubling_solve</a>                                                      -  solves the linear equation X=A*X*B+C
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\lyapunov_equation">m\solvers\X_equal_A_X_B_plus_C_solvers\lyapunov_equation</a>                                                   -  solves the equation V=T*V*T'+Q
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_a_la_tadonki">m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_a_la_tadonki</a>                                               -  attempts to solve the equation V=A*V*B+C
%   <a href="matlab:help m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_solve">m\solvers\X_equal_A_X_B_plus_C_solvers\sandwich_solve</a>                                                      -  solves the linear equation X=A*X*B+C
%   m\solvers\X_equal_A_X_B_plus_C_solvers\tests                                                               - (No help available)
%
%
