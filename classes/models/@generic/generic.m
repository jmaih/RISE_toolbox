classdef generic
    
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
        
        % structure holding information on modifiable settings
        options
        
        % information on estimation: posterior maximization and simulation
        estimation =orderfields(...
        struct('endogenous_priors',[],'priors',[],...
        'posterior_maximization',struct(...
        'estim_start_time',[],'estim_end_time',[],...
        'log_lik',[],'log_post',[],'log_prior',[],'log_endog_prior',[],...
        'active_inequalities_number',0,'hessian',[],'vcov',[],'mode',[],...
        'mode_stdev',[],'funevals',[],...
        'log_marginal_data_density_laplace',[]...
        ),...
        'posterior_simulation',[])...
        );

        % model solution including steady state, definitions, etc.
        solution
        
        % structure holding predicted, updated and smoothed series
        filtering
        
    end
    
    properties(SetAccess=protected)
        % information on markov chains, regimes and related items
        markov_chains
        
    end
        
    properties(SetAccess = protected, Hidden = true)%(SetAccess=protected)
        
        data
        
        data_are_loaded = false

		estim_priors_data = []
        
		estim_endogenous_priors_data = []
        
        estimation_restrictions
        
        estimation_under_way=false
        
        general_restrictions_data
        
        linear_restrictions_data
        
        nonlinear_restrictions_data
        
        list_of_issues
        
        miscellaneous=struct() % will hold elements that are not classified
        
        number_of_restrictions
        
        parameter_values
        
        restrictions_are_absorbed=false
        
        spec_checker
        
    end
    
    properties(Hidden = true)
        
        routines=struct()
        
        online_routines=struct()
        
        % helps return the likelihood instead of the posterior when the
        % computation of the Fisher information matrix is needed.
        is_fisher=false
        
    end
    
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        varargout=conclude_estimation(varargin)
        
    end
    
    methods(Abstract)
        % methods that must be implemented by the subclasses
        varargout=solve(varargin)
        
        % methods that must be implemented by the subclasses
        varargout=set_solution_to_companion(varargin)
        
    end
    
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        
        varargout=problem_reduction(varargin)
        
    end
    
    methods
        
        function obj=generic(varargin)
            
            if nargin
                
                obj=generic.reset(obj,varargin{:});
                
            end
            
        end       
        
        varargout=draw_parameter(varargin)
        
        varargout=evaluate_general_restrictions(varargin)
        
        varargout=estimate(varargin)
        
        varargout=fisher(varargin)
        
        varargout=forecast(varargin)
        
        varargout=get(varargin)
        
        varargout=growth_database(varargin)
        
        varargout=hessian(varargin)
        
        varargout=historical_decomposition(varargin)
        
        varargout=hd(varargin)
        
        varargout=initial_conditions(varargin)
        
        varargout=irf(varargin)
        
        varargout=is_stable_system(varargin)
        
        varargout=isnan(varargin)
        
        varargout=load_parameters(varargin)
        
        varargout=log_posterior_kernel(varargin)
        
        varargout=log_prior_density(varargin)
        
        varargout=mode_curvature(varargin)
        
        varargout=plot_priors(varargin)
        
        varargout=plot_priors_and_posteriors(varargin)
        
        varargout=plot_posteriors(varargin)
        
        varargout=posterior_sample(varargin)
        
        varargout=print_estimation_results(varargin)
        
        varargout=print_estimation_results_legacy(varargin)
        
        varargout=pull_objective(varargin)
        
        varargout=randsample(varargin)
        
        varargout=refresh(varargin)
        
        varargout=report(varargin)
        
        varargout=set(varargin)
        
        varargout=simulate(varargin)
        
        varargout=stoch_simul(varargin)
        
        varargout=theoretical_autocorrelations(varargin)
        
        varargout=theoretical_autocovariances(varargin)
        
        varargout=variance_decomposition(varargin)
        
    end
    
    methods(Static,Hidden=true)
        
        function obj=reset(obj,varargin)
            
            implemented_classes={'svar','rfvar','dsge','vstar'};%,'stochvol'
            
            model_type=class(obj);
            
            endogen=varargin{1};
            
            exogen=varargin{2};
            
            observs=varargin{3};
            
            if ~ismember(model_type,implemented_classes)
                %  var vecm varx favar dfm dsge msdsge msre svar msvar
                %  mssvar dsge-var dsge-kulish star stvar stochvol
                error(['"',model_type,'" model class not yet implemented'])
                
            end
            
            % further declarations
            %----------------------
            obj=do_names(obj,endogen,'endogenous');
            
            obj=do_names(obj,exogen,'exogenous');
            
            obj=do_names(obj,observs,'observables');
            
            mark_parameters()
            
            markchains=varargin{4};
            
            % the markov chains will set the parameters
            %------------------------------------------
            obj=add_markov_chains_and_parameters(obj,markchains);
            
            function mark_parameters()
                
                n=obj.parameters.number;
                
                obj.parameters.is_switching=false(1,n);
                
                obj.parameters.is_trans_prob=false(1,n);
                
                obj.parameters.governing_chain=ones(1,n);
                
                % re-tag the transition probabilities
                %------------------------------------
                for iparam=1:n
                    
                    obj.parameters.is_trans_prob(iparam)=...
                        parser.is_transition_probability(obj.parameters.name{iparam});
                    
                end
                
            end
            
        end
        
        varargout=describe_regimes(varargin)
        
    end
    
    methods(Access=private)
        
        varargout=load_mode(varargin)
        
        varargout=add_markov_chains_and_parameters(varargin)
        
    end
    
    methods(Access=protected)
        
        varargout=parameters_links(varargin)
        
        varargout=setup_priors(varargin)
        
        varargout=setup_endogenous_priors(varargin)
        
        varargout=setup_calibration(varargin)
        
        varargout=setup_general_restrictions(varargin)
        
        varargout=setup_linear_restrictions(varargin)
        
        varargout=setup_nonlinear_restrictions(varargin)
        
    end
    
    methods(Access=protected,Hidden=true)
        
        varargout=add_to_routines(varargin)
        
        varargout=check_property(varargin)
        
        varargout=complementarity_memoizer(varargin)
        
        varargout=do_names(varargin)
        
        varargout=estimation_wrapper(varargin)
        
        varargout=find_posterior_mode(varargin)
        
        varargout=get_estimated_parameter_names(varargin)
        
        varargout=update_estimated_parameter_names(varargin)
        
        varargout=re_order_output_rows(varargin)
        
        varargout=set_simulation_initial_conditions(varargin)
        
        varargout=transform_parameters(varargin)
        
        varargout=unstransform_parameters(varargin)
        
    end
    
    methods(Static,Hidden=true)
        
        varargout=decompose_parameter_name(varargin)
        
    end
    
    methods(Hidden=true)
        
        varargout=prepare_transition_routine(varargin)
        
    end
    
    methods(Hidden=true)
        
        varargout=assign_estimates(varargin)
        
        varargout=back_door(varargin)
        
        varargout=load_data(varargin)
        
        varargout=data_prerequest(varargin)
        
        varargout=quick_irfs(varargin)
        
        varargout=quick_plots(varargin)
        
        varargout=setup_restrictions(varargin)
        
        varargout=stationary_index(varargin)
        
        varargout=unclassified(varargin)
        
    end
    
end

