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
        
		estim_dirichlet=[]
        
        estimation_restrictions
        
        estimation_under_way=false
        
        general_restrictions_data
        
        linear_restrictions_data
        
        nonlinear_restrictions_data
        
        list_of_issues
        
        miscellaneous=struct() % will hold elements that are not classified
        
        number_of_restrictions
        
        parameter_values
        
        routines=struct()
        
        online_routines=struct()
        
        restrictions_are_absorbed=false
        
    end
    
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        varargout=conclude_estimation(varargin)
        
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
        
        varargout=forecast(varargin)
        
        varargout=get(varargin)
        
        varargout=growth_database(varargin)
        
        varargout=hessian(varargin)
        
        varargout=initial_conditions(varargin)
        
        varargout=irf(varargin)
        
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
        
        varargout=pull_objective(varargin)
        
        varargout=randsample(varargin)
        
        varargout=refresh(varargin)
        
        varargout=report(varargin)
        
        varargout=set(varargin)
        
        varargout=simulate(varargin)
        
        varargout=stoch_simul(varargin)
        
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
            
        end
        
    end
    
    methods(Access=private)
        
        varargout=load_mode(varargin)
        
    end
    
    methods(Access=protected)
        
        varargout=decompose_parameter_name(varargin)
        
        varargout=parameters_links(varargin)
        
        varargout=setup_priors(varargin)
        
        varargout=setup_calibration(varargin)
        
        varargout=setup_general_restrictions(varargin)
        
        varargout=setup_linear_restrictions(varargin)
        
        varargout=setup_nonlinear_restrictions(varargin)
        
    end
    
    methods(Access=protected,Hidden=true)
        
        varargout=add_to_routines(varargin)
        
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
    
    methods(Hidden=true)
        
        varargout=assign_estimates(varargin)
        
        varargout=load_data(varargin)
        
        varargout=quick_irfs(varargin)
        
        varargout=quick_plots(varargin)
        
        varargout=setup_restrictions(varargin)
        
        varargout=stationary_index(varargin)
        
    end
    
end

