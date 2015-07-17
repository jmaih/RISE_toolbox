classdef rise_generic % < matlab.mixin.Heterogeneous
    % rise_generic generic class for various types of models
    %
    % methods
    % --------
    %
    % - [check_optimum](rise_generic/check_optimum)
    % - [draw_parameter](rise_generic/draw_parameter)
    % - [estimate](rise_generic/estimate)
    % - [forecast](rise_generic/forecast)
    % - [get](rise_generic/get)
    % - [historical_decomposition](rise_generic/historical_decomposition)
    % - [irf](rise_generic/irf)
    % - [isnan](rise_generic/isnan)
    % - [load_parameters](rise_generic/load_parameters)
    % - [log_marginal_data_density](rise_generic/log_marginal_data_density)
    % - [log_posterior_kernel](rise_generic/log_posterior_kernel)
    % - [log_prior_density](rise_generic/log_prior_density)
    % - [posterior_marginal_and_prior_densities](rise_generic/posterior_marginal_and_prior_densities)
    % - [posterior_simulator](rise_generic/posterior_simulator)
    % - [print_estimation_results](rise_generic/print_estimation_results)
    % - [prior_plots](rise_generic/prior_plots)
    % - [report](rise_generic/report)
    % - [rise_generic](rise_generic/rise_generic)
    % - [set](rise_generic/set)
    % - [set_solution_to_companion](rise_generic/set_solution_to_companion)
    % - [simulate](rise_generic/simulate)
    % - [simulation_diagnostics](rise_generic/simulation_diagnostics)
    % - [solve](rise_generic/solve)
    % - [stoch_simul](rise_generic/stoch_simul)
    % - [theoretical_autocorrelations](rise_generic/theoretical_autocorrelations)
    % - [theoretical_autocovariances](rise_generic/theoretical_autocovariances)
    % - [variance_decomposition](rise_generic/variance_decomposition)
    %
    % properties
    % -----------
    %
    % - [legend] -
    % - [endogenous] -
    % - [exogenous] -
    % - [parameters] -
    % - [observables] -
    % - [markov_chains] -
    % - [options] -
    % - [estimation] -
    % - [solution] -
    % - [filtering] -
    properties
        % attribute for giving a tag to a specific version of a model
        legend='';
    end
    properties(SetAccess=protected)
        % information on endogenous variables (names, number, types, etc.)
        endogenous
        % information on exogenous variables (names, number, types, etc.)
        exogenous
        % information on parameters (names, number, types, etc.)
        parameters
        % information on observable variables (names, number, types, etc.)
        observables
        % information on markov chains, regimes and related items
        markov_chains
        % structure holding information on modifiable settings
        options
        % information on estimation: posterior maximization and simulation
        estimation
        % model solution including steady state, definitions, etc.
        solution
        % structure holding predicted, updated and smoothed series
        filtering
    end
    properties(SetAccess = protected, Hidden = true)%(SetAccess=protected)
        data
        data_are_loaded = false
        estim_distrib_locations={}
        estim_distributions={}
        estim_hyperparams=[]
        estimation_restrictions
        estimation_under_way=false;
        general_restrictions_data
        linear_restrictions_data
        list_of_issues
        miscellaneous=struct() % will hold elements that are not classified
        parameter_values
        routines=struct()
    end
    
    methods(Abstract)
        % methods that must be implemented by the subclasses
        varargout=solve(varargin)
        varargout=set_solution_to_companion(varargin)
    end
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        varargout=conclude_estimation(varargin)
    end
    methods
        function obj=rise_generic(varargin)
            if nargin
                obj=rise_generic.reset(obj,varargin{:});
            end
        end
        varargout=confidence_regions(varargin)
        varargout=draw_parameter(varargin)
        varargout=evaluate_nonlinear_restrictions(varargin)
        varargout=estimate(varargin)
        varargout=forecast(varargin)
        varargout=get(varargin)
        varargout=hessian(varargin)
        varargout=historical_decomposition(varargin)
        varargout=irf(varargin)
        varargout=isnan(varargin)
        varargout=load_parameters(varargin)
        varargout=log_marginal_data_density(varargin)
        varargout=log_posterior_kernel(varargin)
        varargout=log_prior_density(varargin)
        varargout=mode_curvature(varargin)
        varargout=plot_priors(varargin)
        varargout=plot_priors_and_posteriors(varargin)
        varargout=plot_posteriors(varargin)
        varargout=posterior_simulator(varargin)
        varargout=print_estimation_results(varargin)
        varargout=pull_objective(varargin)
        varargout=refresh(varargin)
        varargout=report(varargin)
        varargout=set(varargin)
        varargout=simulation_diagnostics(varargin)
        varargout=simulate(varargin)
        varargout=stoch_simul(varargin)
        varargout=theoretical_autocorrelations(varargin)
        varargout=theoretical_autocovariances(varargin)
        varargout=variance_decomposition(varargin)
    end
    methods(Static,Hidden=true)
        function obj=reset(obj,varargin)
            implemented_classes={'svar','rfvar','dsge','stochvol'};
            model_type=class(obj);
            endogen=varargin{1};
            exogen=varargin{2};
            observs=varargin{3};
            markchains=varargin{4};
            if ~ismember(model_type,implemented_classes)
                %  var vecm varx favar dfm dsge msdsge msre svar msvar
                %  mssvar dsge-var dsge-kulish star stvar stochvol
                error(['"',model_type,'" model class not yet implemented'])
            end
            %             % model type
            %             %-----------
            %             obj.model_class=model_type;
            % further declarations
            %----------------------
            obj=do_names(obj,endogen,'endogenous');
            obj=do_names(obj,exogen,'exogenous');
            obj=do_names(obj,observs,'observables');
            % the markov chains will set the parameters
            %------------------------------------------
            obj=add_markov_chains_and_parameters(obj,markchains);
        end
    end
    methods(Access=private)
        varargout=add_markov_chains_and_parameters(varargin)
        varargout=load_mode(varargin)
    end
    methods(Access=protected)
        varargout=decompose_parameter_name(varargin)
        varargout=parameters_links(varargin)
        varargout=setup_priors(varargin)
        varargout=setup_calibration(varargin)
        varargout=setup_linear_restrictions(varargin)
        varargout=setup_general_restrictions(varargin)
    end
    methods(Access=protected,Hidden=true)
        varargout=complementarity_memoizer(varargin)
        varargout=do_names(varargin)
        varargout=estimation_wrapper(varargin)
        varargout=find_posterior_mode(varargin)
        varargout=get_estimated_parameter_names(varargin)
        varargout=initialize_posterior_simulation(varargin)
        varargout=update_estimated_parameter_names(varargin)
        varargout=update_posterior_simulation_initial_conditions(varargin)
        varargout=prepare_transition_routine(varargin)
        varargout=re_order_output_rows(varargin)
        varargout=set_simulation_initial_conditions(varargin)
    end
    methods(Hidden=true)
        varargout=assign_estimates(varargin)
        varargout=load_data(varargin)
    end
end

