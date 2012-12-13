classdef rise
% Description
% Details
% Usage
% Arguments
% References:

% Author(s)
%   Junior Maih : IMF and dynare team [junior.maih@gmail.com]
% this update:      first version:          
% 
% NOTA BENE:
% 1-All source code available here is distributed without any warranty for
% non-commercial use. 
% 2-There is no implied warranty of Merchantability of fitness for a
% particular purpose. 
% 3-Use at your own risk.
% 4-Should you find something you think is a problem in this code, please
% shoot an email to junior.maih@gmail.com 
    properties
        % initialize this in the proper way in order to avoid problems of
        % concatenation.
        parameters=rise_param.empty(0,0);
    end
    properties (SetAccess = private, Hidden = true)
        current_expansion_order
        current_solution_algorithm
        current_solution_parameters
        data_are_loaded=false;
        dates_filtering % those two options should be moved elsewhere so that they are visible...
        dates_smoothing
        dsge_prior_weight_id
        estimation_restrictions
        estimation_under_way=false;
        estim_hyperparams=[];
        estim_distributions={};
        estim_distrib_locations={};        
        forward_looking_ids
        func_handles
        hybrid_expectations_lambda_id
        is_dsge_var_model
        is_endogenous_switching_model
        is_hybrid_expectations_model
        is_hybrid_model
        is_imposed_steady_state=false;
        is_in_use_parameter
        is_lagrange_multiplier
        is_linear_model
        is_optimal_policy_model
        is_optimal_simple_rule_model
        is_sticky_information_model
        is_purely_backward_looking_model
        is_purely_forward_looking_model
        is_stationary_model
        is_unique_steady_state=false;
        Lead_lag_incidence
        list_of_issues
        measurement_errors_restrictions
        model_derivatives
        orig_endo_names_current
        parameter_restrictions
        parameter_random_draws_restrictions
        planner
        reordering_index
        shadow_transition_matrix
        steady_state_and_balanced_growth_path
        steady_state_file_2_model_communication
        steady_state_shadow_model
        sticky_information_lambda_id
        z_restrictions
    end
    properties(SetAccess=protected)
        % those can be seen but not changed from outside
        options
        definitions=rise_equation.empty(0,0);
        equations=rise_equation.empty(0,0);
        endogenous_priors
        orig_varendo=rise_variable.empty(0,0);
        varendo=rise_variable.empty(0,0);
        varexo=rise_variable.empty(0,0);
        varobs=rise_variable.empty(0,0);
        varobs_exo=rise_variable.empty(0,0);
        cond_varobs=rise_variable.empty(0,0);
        markov_chains
        % structural matrices
        Aplus
        A0
        Aminus
        B
        C
        H % Measurement errors covariance
        W % Loss function weights
        planner_discount
        planner_commitment
        planner_loss
        % solution matrices
        T
        R
        Q % transition matrix
        dsge_var
        % model attributes and properties
        NumberOfEndogenous=[0,0];
        NumberOfExogenous=0;
        NumberOfParameters=0;
        NumberOfObservables=0;
        NumberOfConditionalObservables=[0,0];
        NumberOfRegimes=0;
        NumberOfEquations=0;
        Regimes
        journal='';
        filename='';
        % estimation results
        log_lik
        log_post
        log_prior
        log_endog_prior
        log_mdd_laplace
        log_mdd_mhm
        log_mdd_chib_jeliazkov
        nonlcon_viol_penalty
        estimated_parameters=rise_estim_param.empty(0,0); 
        % I am tempted of moving estimated_parameters to the unprotected
        % properties in order to allow the user to modify the start values
        % from Matlab. But I am afraid a user can misuse this. In any case,
        % the rise_estim_param object is ready for alterations of its
        % start values.
        Hessian
        vcov
% % % % % %         % properties potentially going out
% % % % % %         hist_dec
% % % % % %         counterfactual_parameters
% % % % % %         counterfactual_shocks
        Filters
        state_names
    end
    methods(Hidden=true)
        obj=assign_estimates(obj,params)
    end
    methods
        check_optimum(obj,varargin)
        varargout=counterfactual(obj,varargin)
        varargout=draw_parameter(obj,varargin)
        varargout=estimate(obj,varargin)
        varargout=filter(obj,varargin)
        % the treatment of the two functions below is not satisfactory under multiple
        % regimes. This is something to address
        varargout=evaluate(obj,varargin)
        varargout=forecast(obj,varargin)
        flag=is_stable_system(obj,varargin)
		varargout=load_parameters(obj,varargin)
        varargout=log_marginal_data_density_mhm(obj,varargin)
        varargout=log_marginal_data_density_chib_jeliazkov(obj,varargin)
        varargout=log_posterior_kernel(obj,varargin)
        varargout=log_prior_density(obj,varargin)
        varargout=monte_carlo_filtering(obj,varargin) % <--- formerly sensitivity_analysis
        posterior_marginal_and_prior_densities(obj,varargin)
        varargout=posterior_simulator(obj,varargin)
        print_estimates(obj,varargin)
        print_estimation_results(obj,varargin)
        print_solution(obj,varargin)
        prior_plots(obj,varargin)
        varargout=set_parameters(obj,varargin)%,namesOrIndex,values,regime
        % under estimation, I put the parameters in the second argument
        % (bypass), instead of loading'em from obj.parameters, which given
        % the properties of the rise_param object, will only slow down
        % estimation considerably.
        varargout=simulate(obj,varargin)
        varargout=simulation_diagnostics(obj,varargin)
        varargout=solve(obj,varargin)
        varargout=solve_alternatives(obj,varargin)
        varargout=stoch_simul(obj,varargin)
        varargout=theoretical_autocorrelations(obj,varargin)
        varargout=theoretical_autocovariances(obj,varargin)
        varargout=irf(obj,varargin)
        varargout=set_properties(obj,varargin)
        varargout=set_options(obj,varargin)
        varargout=historical_decomposition(obj,varargin)
        varargout=variance_decomposition(obj,varargin)
        % constructor
        function obj=rise(model_filename,varargin)
            % default options
            if nargin<1
%                 obj=set_options(obj);
                return
            end
            % the constructor needs to have an output
            dictionary=CompileModelFile(model_filename);
            % build the equations object
            quick_fill={'NumberOfEquations','Lead_lag_incidence','filename',...
                'reordering_index','orig_endo_names_current',...
                'planner','is_sticky_information_model',...
                'is_hybrid_expectations_model','is_optimal_policy_model',...
                'is_dsge_var_model','is_optimal_simple_rule_model',...
                'is_hybrid_model','is_purely_backward_looking_model',...
                'is_purely_forward_looking_model','is_lagrange_multiplier',...
                'is_linear_model','shadow_transition_matrix',...
                'is_endogenous_switching_model','is_in_use_parameter'};
            
            % check that all the shocks in the model are in use
            not_in_use={dictionary.exogenous(~dictionary.is_in_use_shock).name};
            if ~isempty(not_in_use)
                disp([mfilename,'(gentle warning) :: the following shocks do not seem to affect the model.'])
                disp('You may want to discard them from your model file for tidiness')
                disp(not_in_use)
            end
            
            for ii=1:numel(quick_fill)
                obj.(quick_fill{ii})=dictionary.(quick_fill{ii});
            end
            obj.steady_state_shadow_model=dictionary.static.steady_state_shadow_model;  
            obj.is_imposed_steady_state=dictionary.static.is_imposed_steady_state;
            obj.is_unique_steady_state=dictionary.static.is_unique_steady_state;
            
            for ii=1:obj.NumberOfEquations
                obj.equations(ii,1)=rise_equation('id',ii,...
                    'dynamic',dictionary.dynamic.model{ii},...
                    'static',dictionary.static.model{ii},...
                    'shadow_dynamic',dictionary.dynamic.shadow_model{ii},...
                    'shadow_static',dictionary.static.shadow_model{ii},...
                    'shadow_balanced_growth_path',char(dictionary.static.shadow_BGP_model(2*(ii-1)+1:2*ii)));
            end
            
            for ii=1:numel(dictionary.definitions)
                obj.definitions(ii,1)=rise_equation('id',ii,...
                    'dynamic',dictionary.definitions(ii).model,...
                    'shadow_dynamic',dictionary.definitions(ii).shadow);
            end
            
            obj.NumberOfEndogenous=[numel(dictionary.orig_endogenous),...
                numel(dictionary.endogenous)];
            obj.NumberOfExogenous=numel(dictionary.exogenous);
            obj.NumberOfParameters=numel(dictionary.parameters);
            obj.NumberOfObservables=numel(dictionary.observables);
            for ii=1:obj.NumberOfEndogenous(1)
                obj.orig_varendo(ii,1)=rise_variable(dictionary.orig_endogenous(ii).name,...
                    'tex_name',dictionary.orig_endogenous(ii).tex_name,'id',ii);
            end
            for ii=1:obj.NumberOfEndogenous(2)
                % HERE WE MAY NEED TO KNOW WHETHER A VARIABLE IS A LAGRANGE
                % MULTIPLIER OR NOT... IS THIS IMPORTANT?
                obj.varendo(ii,1)=rise_variable(dictionary.endogenous(ii).name,...
                    'tex_name',dictionary.endogenous(ii).tex_name,'id',ii);
            end
            for ii=1:obj.NumberOfExogenous
                obj.varexo(ii,1)=rise_variable(dictionary.exogenous(ii).name,...
                    'tex_name',dictionary.exogenous(ii).tex_name,'id',ii);
            end
            
            iter_endo=0;
            iter_exo=0;
            for ii=1:obj.NumberOfObservables
                navn=dictionary.observables(ii).name;
                id=locate_variables(navn,{obj.varendo.name},true);
                % values will be assigned only during estimation. There is
                % of course a question of wether I should compute moments
                % of the observed data. But somehow, this too could be
                % done during estimation, especially if a sub-sample has to
                % be selected.
                if ~isnan(id)
                    iter_endo=iter_endo+1;
                    obj.varobs(iter_endo,1)=rise_variable(navn,...
                        'tex_name',dictionary.endogenous(id).tex_name,'id',id);
                else
                    id=locate_variables(navn,{obj.varexo.name},false);
                    iter_exo=iter_exo+1;
                    obj.varobs_exo(iter_exo,1)=rise_variable(navn,...
                        'tex_name',dictionary.exogenous(id).tex_name,'id',id);
                end
            end
            % now make the number of observables into a vector
            obj.NumberOfObservables=[iter_endo,iter_exo];
            
            % Once the names of the exogenous variables are known, we can
            % build the options. We call the overloaded method
            obj=set_options(obj,varargin{:});
            obj.current_solution_algorithm=obj.options.solver;
            obj.current_expansion_order=obj.options.order;
            
            if isfield(dictionary.dynamic,'model_derivatives')
                obj.model_derivatives=dictionary.dynamic.model_derivatives;
            end
            
            % restrictions for endogenous, shocks and parameters for
            % automatic and numerical differentiation
            ny=nnz(obj.Lead_lag_incidence);
            nx=obj.NumberOfExogenous;
            switching=find([obj.parameters.is_switching]);
            np=numel(switching);
            obj.z_restrictions={1:ny,...
                ny+(1:nx),...
                {switching,ny+nx+(1:np)}};
            
            % parameters and estimated parameters and more ...
            obj=format_parameters(obj,dictionary.parameters,dictionary.Parameterization_block,...
                dictionary.MarkovChains,dictionary.Param_rest_block);
            obj.NumberOfRegimes=size(obj.Regimes,1);
            
            % Now we can load the functions and create the handles
            obj=load_functions(obj);
            
            % flags for the different types of models
            if obj.is_hybrid_expectations_model
                obj.hybrid_expectations_lambda_id=find(strcmp('hybrid_expectations_lambda',{obj.parameters.name}),1);
            end
            
            if obj.is_sticky_information_model
                obj.forward_looking_ids=dictionary.forward_looking_ids;
                obj.sticky_information_lambda_id=find(strcmp('sticky_information_lambda',{obj.parameters.name}),1);
            end
            
            % measurement errors restrictions
            obj.measurement_errors_restrictions=[];
            for ii=1:obj.NumberOfObservables(1) % pick only the endogenous observables
                vi=obj.varobs(ii).name;
                loc=find(strcmp(['stderr_',vi],{obj.parameters.name}));
                % the line above is a bit inefficient. I have to change it
                if ~isempty(loc)
                    obj.measurement_errors_restrictions=...
                        [obj.measurement_errors_restrictions;ii,loc];
                end
            end
            
            % conclude
            if obj.is_sticky_information_model
                disp([mfilename,':: Sticky Information (SI) model detected'])
            elseif obj.is_hybrid_expectations_model
                disp([mfilename,':: Hybrid Expectations (HE) model detected'])
            elseif obj.is_optimal_policy_model
                disp([mfilename,':: Optimal Policy model detected'])
            elseif obj.is_dsge_var_model
                disp([mfilename,':: DSGE-VAR model detected'])
            elseif obj.is_optimal_simple_rule_model
                disp([mfilename,':: Optimal Simple Rule (OSR) model detected'])
            elseif obj.is_endogenous_switching_model
                disp([mfilename,':: Endogenous Switching DSGE model in ',int2str(obj.NumberOfRegimes),' regimes detected'])
            else
                disp([mfilename,':: Exogenous Switching DSGE model in ',int2str(obj.NumberOfRegimes),' regimes detected'])
            end
            if obj.is_hybrid_model
                disp([mfilename,':: model has both backward and forward-looking components'])
            elseif obj.is_purely_backward_looking_model
                disp([mfilename,':: model is purely backward looking'])
            elseif obj.is_purely_forward_looking_model
                disp([mfilename,':: model is purely forward-looking'])
            else
                disp([mfilename,':: model is purely static'])
            end
        end
    end
end

