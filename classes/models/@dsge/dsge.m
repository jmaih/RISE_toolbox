classdef dsge < rise_generic
    % dsge Markov Switching Dynamic Stochastic General Equilibrium Modeling
    %
    % dsge Methods:
    % ----------------
    %
    % check_derivatives -  compares the derivatives and the solutions from various differentiation techniques
    % check_optimum -   H1 line
    % compute_steady_state -   H1 line
    % create_estimation_blocks -   H1 line
    % create_state_list - creates the list of the state variables in the solution
    % draw_parameter -   H1 line
    % dsge -   default options
    % estimate -  estimates the parameters of a RISE model
    % filter -   H1 line
    % forecast -  computes forecasts for rise|dsge|svar|rfvar models
    % forecast_real_time -  forecast from each point in time
    % get -   H1 line
    % historical_decomposition - Computes historical decompositions of a DSGE model
    % irf -   H1 line
    % is_stable_system -   H1 line
    % isnan -   H1 line
    % load_parameters -   H1 line
    % log_marginal_data_density -   H1 line
    % log_posterior_kernel -   H1 line
    % log_prior_density -   H1 line
    % monte_carlo_filtering -   H1 line
    % posterior_marginal_and_prior_densities -   H1 line
    % posterior_simulator -   H1 line
    % print_estimation_results -   H1 line
    % print_solution -  print the solution of a model or vector of models
    % prior_plots -   H1 line
    % refresh -  refresh the options of an old object with a newer version of
    % report - assigns the elements of interest to a rise_report.report object
    % resid -   H1 line
    % set -  sets options for dsge|rise models
    % set_solution_to_companion -   H1 line
    % simulate -  simulates a RISE model
    % simulate_nonlinear -   H1 line
    % simulation_diagnostics -   H1 line
    % solve -   H1 line
    % solve_alternatives -   H1 line
    % stoch_simul -   H1 line
    % theoretical_autocorrelations -   H1 line
    % theoretical_autocovariances -   H1 line
    % variance_decomposition -   H1 line
    %
    % dsge Properties:
    % -------------------
    %
    % definitions -   values of auxiliary parameters defined in the model file with a #
    % equations - of the system
    % folders_paths -   paths for the different folders in which RISE stores information
    % dsge_var -
    % filename -   name of the rs/rz/dsge file read
    % legend -   attribute for giving a tag to a specific version of a model
    % endogenous -   information on endogenous variables (names, number, types, etc.)
    % exogenous -   information on exogenous variables (names, number, types, etc.)
    % parameters -   information on parameters (names, number, types, etc.)
    % observables -   information on observable variables (names, number, types, etc.)
    % markov_chains -   information on markov chains, regimes and related items
    % options -   structure holding information on modifiable settings
    % estimation -   information on estimation: posterior maximization and simulation
    % solution -   model solution including steady state, definitions, etc.
    % filtering -   structure holding predicted, updated and smoothed series
    properties (Hidden = true)
        online_routines
    end
    properties (SetAccess = private, Hidden = true)
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
        switching_parameters_leads_index
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
        % provision for functions, solving and resolving
        %--------------------------------------------------
        warrant_resolving = true;
        warrant_setup_change = true % initialization of functions, derivatives, etc.
    end
    properties(SetAccess=protected)
        % values of auxiliary parameters defined in the model file with a #
        definitions
        % equations of the system
        equations
        % paths for the different folders in which RISE stores information
        folders_paths
        dsge_var
        % name of the rs/rz/dsge file read
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
            % make a copy of the routines
            %-----------------------------
            obj.online_routines=obj.routines;
            
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
            
            % Then load the functions/routines: some routines may be
            % created inside format_parameters
            %--------------------------------------------------------------
            create_folders_and_add_further_routines();
            
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
                    likelihood=@likelihood_dsge_var;
                elseif obj.is_optimal_simple_rule_model
                    likelihood=@likelihood_optimal_simple_rule;
                else
                    likelihood=@likelihood_markov_switching_dsge;
                end
                obj=add_to_routines(obj,'likelihood',likelihood);
            end
        end
    end
    methods(Sealed)
        varargout=simulate(varargin)
    end
    methods(Access=protected)
        varargout=setup_calibration(varargin)
    end
    methods(Access=protected,Hidden=true)
        varargout=re_order_output_rows(varargin)
        varargout=prepare_transition_routine(varargin)
    end
    methods(Hidden=true)
        varargout=assign_estimates(varargin)
        varargout=do_not_anticipate_future_shocks(varargin)
        varargout=filter_initialization(varargin)
        varargout=latex_model_file(varargin)
        varargout=load_solution(varargin)
        varargout=load_data(varargin)
        varargout=set_z_eplus_horizon(varargin)
    end
end

