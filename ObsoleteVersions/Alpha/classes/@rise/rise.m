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
        legend='' 
        % could be useful for differentiating models in reports. readily
        % accessible through obj.legend
    end
    properties (SetAccess = private, Hidden = true)
        current_expansion_order
        current_solution_algorithm
        current_solution_parameters
        current_solve_order
        data_are_loaded=false;
        dates_filtering % those two options should be moved elsewhere so that they are visible...
        dates_smoothing
        dsge_prior_weight_id
        estimation_restrictions
        estimation_under_way=false;
        estim_hyperparams=[];
        estim_distributions={};
        estim_distrib_locations={};
        % structure
        evaluated_derivatives=struct('zeroth_order',{},'first_order',{},'second_order',{})
        forward_looking_ids
        func_handles
%         G0
        G1 % holder of packed first-order derivatives
        G2 % holder of packed second-order derivatives
        hybrid_expectations_lambda_id
        is_dsge_var_model
        is_endogenous_switching_model
        is_hybrid_expectations_model
        is_hybrid_model
        is_imposed_steady_state=false;
        is_in_use_parameter
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
        data=struct() % structure to hold data, solution matrices for estimation
        measurement_errors_restrictions
        model_derivatives
        orig_endo_names_current
        parameter_values=[] % nparam x nregimes matrix
        parameter_restrictions
        planner_system
        raw_file
        rawfile_triggers
        reordering_index
        shadow_transition_matrix
        steady_state_file_2_model_communication
        steady_state_shadow_model
        sticky_information_lambda_id
        input_list
        z_restrictions
    end
    properties(SetAccess=protected)
        % those can be seen but not changed from outside at least not directly
        % initialize this in the proper way in order to avoid problems of
        % concatenation.
        options
        definitions
        equations
        folders_paths
        markov_chains
        endogenous=struct()
        exogenous=struct()
        parameters=struct()
        observables=struct()
% %         planner=struct('discount',{},'commitment',{},'objective',{},'weights',{})
% %         % first and second-order approximation
% %         structural
        % solution matrices
        solution
        dsge_var
        % model attributes and properties
        filename='';
        filtering
        estimation
    end
    methods(Hidden=true)
        varargout=assign_estimates(obj,varargin)
    end
    methods
        varargout=check_optimum(obj,varargin)
        varargout=counterfactual(obj,varargin)
        varargout=check_derivatives(obj,varargin)
        varargout=draw_parameter(obj,varargin)
        varargout=estimate(obj,varargin)
        varargout=filter(obj,varargin)
        % the treatment of the two functions below is not satisfactory under multiple
        % regimes. This is something to address
        varargout=forecast(obj,varargin)
        varargout=forecast_real_time(obj,varargin)
        varargout=get(obj,varargin)
        varargout=isnan(varargin)
        varargout=is_stable_system(obj,varargin)
        varargout=load_parameters(obj,varargin)
        varargout=log_marginal_data_density(obj,varargin)
        varargout=log_posterior_kernel(obj,varargin)
        varargout=log_prior_density(obj,varargin)
        varargout=monte_carlo_filtering(obj,varargin) % <--- formerly sensitivity_analysis
        varargout=posterior_marginal_and_prior_densities(obj,varargin)
        varargout=posterior_simulator(obj,varargin)
        varargout=print_estimation_results(obj,varargin)
        varargout=print_solution(obj,varargin)
        varargout=prior_plots(varargin)
        varargout=report(obj,varargin)
        varargout=set(obj,varargin)
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
                return
            end
            if isa(model_filename,'rise')
                obj=model_filename;
                return
            end
            
            % separate options going into CompileModelFile to the others
            cmfOptions=fieldnames(parser.parse());
            if rem(length(varargin),2)~=0
                error('arguments must come in pairs')
            end
            discard=false(1,length(varargin));
            for iv=1:2:length(varargin)
                if ~ischar(varargin{iv})
                    error('odd input arguments must be char')
                end
                if ismember(varargin{iv},cmfOptions)
                    discard(iv:iv+1)=true;
                end
            end
            cmfArgins=varargin(discard);
            % proceed with the remaining arguments
            varargin(discard)=[];
            % the constructor needs to have an output
            dictionary=parser.parse(model_filename,cmfArgins{:});
            % build the equations object
            quick_fill={'Lead_lag_incidence','filename',...
                'reordering_index','orig_endo_names_current',... %
                'planner_system','is_sticky_information_model',...
                'is_hybrid_expectations_model','is_optimal_policy_model',...
                'is_dsge_var_model','is_optimal_simple_rule_model',...
                'is_hybrid_model','is_purely_backward_looking_model',...
                'is_purely_forward_looking_model',...
                'is_linear_model','shadow_transition_matrix',...
                'is_endogenous_switching_model','input_list',...
                'is_imposed_steady_state','is_unique_steady_state',...
                'endogenous','exogenous','parameters','observables',...
                'raw_file','rawfile_triggers','equations','definitions',...
                'model_derivatives','steady_state_shadow_model',...
                'markov_chains'};
            
            % check that all the shocks in the model are in use
            not_in_use=dictionary.exogenous.name(~dictionary.exogenous.is_in_use);
            if ~isempty(not_in_use)
                disp([mfilename,'(gentle warning) :: the following shocks do not seem to affect the model.'])
                disp('You may want to discard them from your model file for tidiness')
                disp(not_in_use)
            end
            
            for ii=1:numel(quick_fill)
                obj.(quick_fill{ii})=dictionary.(quick_fill{ii});
            end                        
            
            % flags for the different types of models
            if obj.is_hybrid_expectations_model
                obj.hybrid_expectations_lambda_id=find(strcmp('hybrid_expectations_lambda',{obj.parameters.name}),1);
            end
            
            if obj.is_sticky_information_model
                obj.forward_looking_ids=dictionary.forward_looking_ids;
                obj.sticky_information_lambda_id=find(strcmp('sticky_information_lambda',{obj.parameters.name}),1);
            end
            
            % Once the names of the exogenous variables are known, we can
            % build the options. We call the overloaded method
            obj=set_options(obj,varargin{:});
            % which will also load the some of the functions since
            % obj.options is empty at this point.
            
            % parameters and estimated parameters and more ...:
            % format_parameters depends on the value of the prior
            % truncation which is stored in the options... and cannot be
            % changed after the model is read. yet another reason to
            % separate the reading of the model with everything else
            obj=format_parameters(obj,dictionary.Parameterization_block,...
                dictionary.Param_rest_block);
            
            % initialize output matrices
%             obj.structural=rise.initialize_solution_or_structure('system',obj.markov_chains.regimes_number);
            obj.solution=rise.initialize_solution_or_structure('solution',obj.markov_chains.regimes_number);
            
            % restrictions for endogenous, shocks and parameters for
            % automatic and numerical differentiation
            ny=nnz(obj.Lead_lag_incidence);
            nx=sum(obj.exogenous.number);
            switching=obj.parameters.is_switching;
            np=sum(switching);
            obj.z_restrictions={1:ny,...
                ny+(1:nx),...
                {switching,ny+nx+(1:np)}};
            
            % the two lines below should probably be moved elsewhere
            obj.current_solution_algorithm=obj.options.solver;
            obj.current_expansion_order=obj.options.solve_expect_order;
            
            % conclude
            disp(' ')
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
                disp([mfilename,':: Endogenous Switching DSGE model in ',int2str(obj.markov_chains.regimes_number),' regimes detected'])
            else
                disp([mfilename,':: Exogenous Switching DSGE model in ',int2str(obj.markov_chains.regimes_number),' regimes detected'])
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
    methods(Static,Hidden=true)
        varargout=initialize_solution_or_structure(varargin)
    end
    methods(Static)
    end
end

