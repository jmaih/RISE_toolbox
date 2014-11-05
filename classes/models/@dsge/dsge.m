classdef dsge < rise_generic
    % dsge Markov Switching Dynamic Stochastic General Equilibrium Modeling
    %
    % methods
    % --------
    %
    % - [check_derivatives](dsge/check_derivatives)
    % - [check_optimum](dsge/check_optimum)
    % - [compute_steady_state](dsge/compute_steady_state)
    % - [create_estimation_blocks](dsge/create_estimation_blocks)
    % - [create_state_list](dsge/create_state_list)
    % - [draw_parameter](dsge/draw_parameter)
    % - [dsge](dsge/dsge)
    % - [estimate](dsge/estimate)
    % - [filter](dsge/filter)
    % - [forecast](dsge/forecast)
    % - [forecast_real_time](dsge/forecast_real_time)
    % - [get](dsge/get)
    % - [historical_decomposition](dsge/historical_decomposition)
    % - [irf](dsge/irf)
    % - [is_stable_system](dsge/is_stable_system)
    % - [isnan](dsge/isnan)
    % - [load_parameters](dsge/load_parameters)
    % - [log_marginal_data_density](dsge/log_marginal_data_density)
    % - [log_posterior_kernel](dsge/log_posterior_kernel)
    % - [log_prior_density](dsge/log_prior_density)
    % - [monte_carlo_filtering](dsge/monte_carlo_filtering)
    % - [posterior_marginal_and_prior_densities](dsge/posterior_marginal_and_prior_densities)
    % - [posterior_simulator](dsge/posterior_simulator)
    % - [print_estimation_results](dsge/print_estimation_results)
    % - [print_solution](dsge/print_solution)
    % - [prior_plots](dsge/prior_plots)
    % - [report](dsge/report)
    % - [resid](dsge/resid)
    % - [set](dsge/set)
    % - [set_solution_to_companion](dsge/set_solution_to_companion)
    % - [simulate](dsge/simulate)
    % - [simulate_nonlinear](dsge/simulate_nonlinear)
    % - [simulation_diagnostics](dsge/simulation_diagnostics)
    % - [solve](dsge/solve)
    % - [solve_alternatives](dsge/solve_alternatives)
    % - [stoch_simul](dsge/stoch_simul)
    % - [theoretical_autocorrelations](dsge/theoretical_autocorrelations)
    % - [theoretical_autocovariances](dsge/theoretical_autocovariances)
    % - [variance_decomposition](dsge/variance_decomposition)
    %
    % properties
    % -----------
    %
    % - [definitions] -
    % - [equations] -
    % - [folders_paths] -
    % - [dsge_var] -
    % - [filename] -
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
    properties (Hidden = true)
        online_routines
        disc_routines
    end
    properties (SetAccess = private, Hidden = true)
        current_solution_state
        dates_filtering % those two options should be moved elsewhere so that they are visible...
        dates_smoothing
        dsge_prior_weight_id
        forward_looking_ids
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
        lead_lag_incidence
        measurement_errors_restrictions
        model_derivatives
        number_of_restrictions
        orig_endo_names_current
        planner_system
        raw_file
        rawfile_triggers
        steady_state_file_2_model_communication
        steady_state_funcs % holder for steady-state calculating functions
        sticky_information_lambda_id
        input_list
        % statistics on the order of the variables
        %--------------------------------------------
        v
        locations
        siz
        order_var
        inv_order_var
        % provision for automatic differentiation
        %----------------------------------------
        steady_state_index
        % provision for loose commitment and sticky information
        %------------------------------------------------------
        reordering_index
        equations_reordering_for_multipliers
    end
    properties(SetAccess=protected)
        % those can be seen but not changed from outside at least not directly
        % initialize this in the proper way in order to avoid problems of
        % concatenation.
        definitions
        equations
        folders_paths
        dsge_var
        filename='';
    end
    methods(Hidden)
        varargout=conclude_estimation(varargin)
    end
    methods
        varargout=check_derivatives(varargin)
        varargout=create_state_list(varargin)
        varargout=filter(varargin)
        varargout=print_solution(varargin)
        varargout=solve(varargin)
        varargout=solve_alternatives(varargin)
        % the treatment of the two functions below is not satisfactory under multiple
        % regimes. This is something to address
        varargout=forecast_real_time(varargin)
        varargout=is_stable_system(varargin)
        varargout=monte_carlo_filtering(varargin)
        varargout=resid(varargin)
        varargout=simulate_nonlinear(varargin)
        varargout=set_solution_to_companion(varargin)
        % constructor
        %------------
        function obj=dsge(model_filename,varargin)
            % default options
            obj=obj@rise_generic();
            if nargin<1
                return
            elseif isa(model_filename,'dsge')
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
            quick_fill={'lead_lag_incidence','filename',...
                'reordering_index','orig_endo_names_current',... %
                'planner_system','is_sticky_information_model',...
                'is_hybrid_expectations_model','is_optimal_policy_model',...
                'is_dsge_var_model','is_optimal_simple_rule_model',...
                'is_hybrid_model','is_purely_backward_looking_model',...
                'is_purely_forward_looking_model',...
                'is_linear_model',...
                'is_endogenous_switching_model','input_list',...
                'is_imposed_steady_state','is_unique_steady_state',...
                'endogenous','exogenous','parameters','observables',...
                'raw_file','rawfile_triggers','equations','definitions',...
                'markov_chains','v','locations','siz','order_var',...
                'inv_order_var','steady_state_index',...
                'equations_reordering_for_multipliers','routines'};
            
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
                obj.hybrid_expectations_lambda_id=find(strcmp('hybrid_expectations_lambda',obj.parameters.name),1);
            end
            
            if obj.is_sticky_information_model
                obj.forward_looking_ids=dictionary.forward_looking_ids;
                obj.sticky_information_lambda_id=find(strcmp('sticky_information_lambda',obj.parameters.name),1);
            end
            
            % Once the names of the exogenous variables are known, we can
            % build the options.
            %--------------------------------------------------------------
            obj=set(obj,varargin{:});
            % Then load the functions/routines: should consider doing this
            % in the parser. The question is whether there is anything that
            % might need re-updating in case options change, like for
            % instance the main folder...
            %--------------------------------------------------------------
            create_folders_and_add_further_routines();
            
            if obj.options.debug
                dicfields=fieldnames(dictionary);
                superfluous=setdiff(dicfields,quick_fill);
                if ~isempty(superfluous)
                    disp(superfluous(:)')
                    warning('the dictionary fields above are superfluous for the constructor')
                end
            end
            
            % parameters and estimated parameters and more ...:
            % format_parameters depends on the value of the prior
            % truncation which is stored in the options... and cannot be
            % changed after the model is read. yet another reason to
            % separate the reading of the model with everything else
            obj=format_parameters(obj,dictionary.Parameterization_block,...
                dictionary.Param_rest_block);
            
            % conclude
            disp(' ')
            if obj.is_sticky_information_model
                disp([mfilename,':: Sticky Information (SI) model detected'])
            elseif obj.is_hybrid_expectations_model
                disp([mfilename,':: Hybrid Expectations (HE) model detected'])
            elseif obj.is_optimal_policy_model
                disp([mfilename,':: Switching Optimal Policy model in ',int2str(obj.markov_chains.regimes_number),' regimes detected'])
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
            
            function create_folders_and_add_further_routines()
                MainFolder=obj.options.results_folder;
                SubFoldersList={'graphs','estimation','simulations','routines'};
                
                if ~exist(MainFolder,'dir')
                    mkdir(MainFolder)
                end
                obj.folders_paths=struct();
                for ifold=1:numel(SubFoldersList)
                    subfolder=[MainFolder,filesep,SubFoldersList{ifold}];
                    obj.folders_paths.(SubFoldersList{ifold})=subfolder;
                    if ~exist(subfolder,'dir')
                        mkdir(subfolder)
                    end
                end
                % likelihood functions
                %---------------------
                if obj.is_dsge_var_model
                    obj.routines.likelihood=@likelihood_dsge_var;
                elseif obj.is_optimal_simple_rule_model
                    obj.routines.likelihood=@likelihood_optimal_simple_rule;
                else
                    obj.routines.likelihood=@likelihood_markov_switching_dsge;
                end
                % initialize the holders for routines in case of swap between online and disc
                obj.online_routines=[];
                obj.disc_routines=[];
            end
        end
    end
    methods(Sealed)
        varargout=simulate(varargin)
    end
    methods(Access=protected,Hidden=true)
        varargout=re_order_output_rows(varargin)
        varargout=prepare_transition_routine(varargin)
    end
    methods(Hidden=true)
        varargout=load_data(varargin)
        varargout=do_not_anticipate_future_shocks(varargin)
        varargout=set_z_eplus_horizon(varargin)
        varargout=latex_model_file(varargin)
        varargout=load_solution(varargin)
    end
end

