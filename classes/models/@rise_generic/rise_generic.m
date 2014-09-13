classdef rise_generic % < matlab.mixin.Heterogeneous
    properties
        legend='';
    end
    properties(SetAccess=protected)
        endogenous
        exogenous
        parameters
        observables
        markov_chains
        options
        estimation
        solution
        filtering
    end
    properties(SetAccess = protected, Hidden = true)%(SetAccess=protected)
        estimation_restrictions
        parameter_values
        data
        data_are_loaded = false
        estimation_under_way=false;
        estim_hyperparams=[];
        estim_distributions={};
        estim_distrib_locations={};
        linear_restrictions_data
        general_restrictions_data
        list_of_issues
        routines=struct();
    end
    properties(SetAccess = private)%(SetAccess=protected), Hidden = true
    end
    methods(Abstract)
        % methods that must be implemented by the subclasses
        varargout=solve(varargin)
        varargout=set_solution_to_companion(varargin)
    end
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        varargout=conclude_estimation(varargin)
        varargout=load_order_var_solution(varargin)
    end
    methods
        function obj=rise_generic(varargin)
            if nargin
                obj=rise_generic.reset(obj,varargin{:});
            end
        end
        varargout=estimate(varargin)
        varargout=irf(varargin)
        varargout=simulate(varargin)
        varargout=forecast(varargin)
        varargout=get(varargin)
        varargout=set(varargin)
        varargout=log_prior_density(varargin)
        varargout=prior_plots(varargin)
        varargout=variance_decomposition(varargin)
        varargout=isnan(varargin)
        varargout=draw_parameter(varargin)
        varargout=check_optimum(varargin)
        varargout=load_parameters(varargin)
        varargout=theoretical_autocorrelations(varargin)
        varargout=theoretical_autocovariances(varargin)
        varargout=historical_decomposition(varargin)
        varargout=report(varargin)
        varargout=simulation_diagnostics(varargin)
        varargout=posterior_marginal_and_prior_densities(varargin)
        varargout=posterior_simulator(varargin)
        varargout=log_posterior_kernel(varargin)
        varargout=log_marginal_data_density(varargin)
        varargout=stoch_simul(varargin)
        varargout=print_estimation_results(varargin)
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
        varargout=do_names(varargin)
        varargout=get_estimated_parameter_names(varargin)
        varargout=update_estimated_parameter_names(varargin)
        varargout=find_posterior_mode(varargin)
        varargout=estimation_wrapper(varargin)
        varargout=update_posterior_simulation_initial_conditions(varargin)
        varargout=initialize_posterior_simulation(varargin)
    end
    methods(Hidden=true)
        varargout=load_data(varargin)
        varargout=assign_estimates(varargin)
    end
end

