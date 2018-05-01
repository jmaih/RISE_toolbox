classdef dsge < generic % & gogetter
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
% posterior_marginal_and_prior_densities -   H1 line
% print_estimation_results -   H1 line
% print_solution -  print the solution of a model or vector of models
% print_solution_legacy -  old form of print_solution
% prior_plots -   H1 line
% refresh -  refresh the options of an old object with a newer version of
% report - assigns the elements of interest to a rise_report.report object
% resid -   H1 line
% set -  sets options for dsge|rise models
% set_solution_to_companion -   H1 line
% simulate -  simulates a RISE model
% simulate_nonlinear -   H1 line
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
properties
% this is so that the user can tailor the information to pass
% around functions written by him but called by RISE
user_data

end

properties (Hidden = true)

old_solution

end

properties (Dependent,Hidden)

log_vars

end

properties (SetAccess = private, Hidden = true)
auxiliary_variables % variables for which the user does not need to solve for the steady state

dates_filtering % those two options should be moved elsewhere so that they are visible...

dates_smoothing

hybrid_expectations_lambda_id

is_dsge_var_model

is_endogenous_switching_model

is_hybrid_expectations_model

is_hybrid_model

is_in_use_parameter

is_optimal_policy_model

is_optimal_simple_rule_model

is_purely_backward_looking_model

is_purely_forward_looking_model

lead_lag_incidence

measurement_errors_restrictions

model_derivatives

planner_system

raw_file

rawfile_triggers

sticky_information_lambda_id

input_list

% steady state solution facilitators
%------------------------------------
occurrence

fast_sstate_occurrence

steady_state_blocks

steady_state_2_model_communication

steady_state_2_blocks_optimization

% steady state model characteristics
%------------------------------------
is_param_changed_in_ssmodel

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

% provision for functions, solving and resolving
%--------------------------------------------------
warrant_resolving = true;

warrant_setup_change = true % initialization of functions, derivatives, etc.

% dsge_var stuff
%---------------
dsge_var

dsge_prior_weight_id

% user information
%------------------
user_endogenous_priors_info

end

properties(SetAccess=protected)
% values of auxiliary parameters defined in the model file with a #
definitions

% equations of the system
equations

% paths for the different folders in which RISE stores information
folders_paths

% name of the rs/rz/dsge file read
filename='';

end

methods
% constructor
%------------
function obj=dsge(model_filename,varargin)
% dsge -- constructor for dsge models
%
% ::
%
%
%   obj=dsge(model_filename,varargin)
%
% Args:
%              %
%              % - **model_filename** [char]: name of the model file. The file
%              % should have extensions rs, rz or dsge
%              %
%              % - **varargin** []: pairwise arguments with possiblities as
%              % follows:
%              %
%              %   - **parameter_differentiation** [true|{false}]: compute or
%              %   not parameter derivatives
%              %
%              %   - **definitions_inserted** [true|{false}]: substitute
%              %   definitions given in the model block. Necessary if the
%              %   definitions contain variables.
%              %
%              %   - **definitions_in_param_differentiation** [true|{false}]:
%              %   insert or not definitions in equations before
%              %   differentiating with respect to parameters
%              %
%              %   - **rise_save_macro** [true|{false}]: save the macro file
%              %   in case of insertions of sub-files
%              %
%              %   - **max_deriv_order** [integer|{2}]: order for symbolic
%              %   differentiation. It is recommended to set to 1, especially
%              %   for large models in case one does not intend to solve
%              %   higher-order approximations
%              %
%              %   - **parse_debug** [true|{false}]: debugging in the parser
%              %
%              %   - **add_welfare** [true|{false}]: add the welfare equation
%              %   when doing optimal policy. N.B: The welfare variable, WELF
%              %   is the true welfare multiplied by (1-discount). The
%              %   within-period utility variable, UTIL is automatically
%              %   added. The reason why welfare is not automatically added
%              %   is that oftentimes models with that equation do not solve.
%              %
%              %   - **rise_flags** [struct|cell]: instructions for the
%              %   partial parsing of the rise file. In case of a cell, the
%              %   cell should be a k x 2 cell, where the first column
%              %   collects the conditional parsing names and the second
%              %   column the values.
%              %
% Returns:
%    :
%              %
%              % - **obj** [rise|dsge]: model object
%              %
% Note:
%              %
%              % - The pairwise options listed above are the ones that will be
%              % processed in the parser. Additional options related to
%              % specific methods can also be passed at this stage, but will
%              % only be applied or used when the specific method dealing with
%              % them is called.
%              %
%              % - In RISE it is possible to declare exogenous and make them
%              % observable at the same time. The exogenous that are observed
%              % are determisitic. This is the way to introduce e.g. time
%              % trends. This strategy also opens the door for estimating
%              % partial equilibrium models.
%              %
% Example:
            %
            % See also:
            
            obj=obj@generic();
            
            % obj=obj@gogetter();
            
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
                'planner_system',...
                'auxiliary_variables',...
                'is_hybrid_expectations_model','is_optimal_policy_model',...
                'is_dsge_var_model','is_optimal_simple_rule_model',...
                'is_hybrid_model','is_purely_backward_looking_model',...
                'is_purely_forward_looking_model',...
                'is_endogenous_switching_model',...
                'input_list','is_param_changed_in_ssmodel',...
                'endogenous','exogenous','parameters','observables',...
                'raw_file','rawfile_triggers','equations','definitions',...
                'markov_chains','v','locations','siz','order_var',...
                'inv_order_var','steady_state_index','occurrence',...
                'fast_sstate_occurrence','routines'};
            
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
            
            obj.options.estim_nonlinear_restrictions=dictionary.Param_rest_block;
            
            % parameters and estimated parameters and more ...:
            % format_parameters depends on the value of the prior
            % truncation which is stored in the options... and cannot be
            % changed after the model is read. yet another reason to
            % separate the reading of the model with everything else
            obj=format_parameters(obj,dictionary.Parameterization_block);
            
            % Then load the functions/routines: some routines may be
            % created inside format_parameters
            %--------------------------------------------------------------
            create_folders_and_add_further_routines();
            
            % conclude
            disp(' ')
            
            switch_type='Exogenous';
            
            if obj.is_endogenous_switching_model
                
                switch_type='Endogenous';
                
            end
            
            if obj.is_optimal_policy_model
                
                model_type='Optimal Policy';
                
            elseif obj.is_dsge_var_model
                
                model_type='DSGE-VAR';
                
            elseif obj.is_optimal_simple_rule_model
                
                model_type='Optimal Simple Rule (OSR)';
                
            else
                
                model_type='DSGE';
                
            end
            
            disp([mfilename,':: ',switch_type,' Switching ',model_type,...
                ' model in ',int2str(obj.markov_chains.regimes_number),...
                ' regimes detected'])
            
            if obj.is_hybrid_model
                
                disp([mfilename,':: model has both backward and forward-looking components'])
                
            elseif obj.is_purely_backward_looking_model
                
                disp([mfilename,':: model is purely backward looking'])
                
            elseif obj.is_purely_forward_looking_model
                
                disp([mfilename,':: model is purely forward-looking'])
                
            else
                
                disp([mfilename,':: model is purely static'])
                
            end
            
            if obj.is_hybrid_expectations_model
                
                disp([mfilename,':: The model features Hybrid Expectations'])
                
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
        
        function w=get.log_vars(obj)
            
            w=obj.endogenous.is_log_var|obj.endogenous.is_log_expanded;
            
        end
        
        varargout=bvar_dsge(varargin)
        
        varargout=check_derivatives(varargin)
        
        varargout=create_state_list(varargin)
        
        varargout=estimate(varargin) % re-signed...
        
        varargout=filter(varargin)
        
        varargout=filter_initialization(varargin)
        
        varargout=forecast_real_time(varargin)
        
        varargout=frontier(varargin)
        
% % % %         varargout=get(varargin): TODO
        varargout=is_stationary_system(varargin)
        
        varargout=loss(varargin)
        
        varargout=observables_decomposition(varargin)
        
        varargout=pull_objective(varargin)
        
        varargout=print_solution(varargin)
        
        varargout=print_solution_legacy(varargin)
        
        varargout=refresh(varargin)
        
        varargout=resid(varargin)
        
        varargout=set_solution_to_companion(varargin)
        
        varargout=simulate_nonlinear(varargin)
        
        varargout=solve(varargin)
        
        varargout=solve_alternatives(varargin)
        
    end
    
    methods(Sealed)
        
        varargout=simulate(varargin)
        
    end
    
    methods(Access=protected)
        
        varargout=setup_calibration(varargin)
        
    end
    
    methods(Access=protected,Hidden=true)
        
        varargout=re_order_output_rows(varargin)
        
        varargout=check_property(varargin)
        
    end
    
    methods(Hidden=true)
        
        varargout=assign_estimates(varargin)
        
        varargout=conclude_estimation(varargin) % abstract method
        
        varargout=latex_model_file(varargin)
        
        varargout=load_solution(varargin)
        
        varargout=load_data(varargin)
        
        varargout=prepare_transition_routine(varargin)
        
        varargout=problem_reduction(varargin)
        
        varargout=riff_erize(varargin)
        
        varargout=set_z_eplus_horizon(varargin)
        
        varargout=growth_component_solver(varargin)
        
        varargout=set_planner_derivatives(varargin)
        
    end
    
    methods(Static,Hidden)
        
        varargout=state_space_wrapper(varargin)
        
    end
    
end

