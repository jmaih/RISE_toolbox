classdef svar < rise_generic
    % svar Structural VAR modeling
    %
    % svar Methods:
    % ----------------
    %
    % check_optimum -   H1 line
    % draw_parameter -   H1 line
    % estimate -  estimates the parameters of a RISE model
    % forecast -  computes forecasts for rise|dsge|svar|rfvar models
    % get -   H1 line
    % historical_decomposition - Computes historical decompositions of a DSGE model
    % irf -  computes impulse responses for a RISE model
    % isnan -   H1 line
    % load_parameters -   H1 line
    % log_marginal_data_density -   H1 line
    % log_posterior_kernel -   H1 line
    % log_prior_density -   H1 line
    % msvar_priors -   H1 line
    % posterior_marginal_and_prior_densities -   H1 line
    % posterior_simulator -   H1 line
    % print_estimation_results -   H1 line
    % prior_plots -   H1 line
    % refresh -  refresh the options of an old object with a newer version of
    % report - assigns the elements of interest to a rise_report.report object
    % set -  sets options for RISE models
    % set_solution_to_companion -   H1 line
    % simulate -  simulates a RISE model
    % simulation_diagnostics -   H1 line
    % solve -   H1 line
    % stoch_simul -   H1 line
    % svar - Structural VAR modeling
    % template -
    % theoretical_autocorrelations -   H1 line
    % theoretical_autocovariances -   H1 line
    % variance_decomposition -   H1 line
    %
    % svar Properties:
    % -------------------
    %
    % constant -   true if VAR has a constant
    % nlags -   number of lags in the VAR
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
    end
    properties(SetAccess=protected)
        % true if VAR has a constant
        constant
        % number of lags in the VAR
        nlags
    end
    properties(SetAccess = private, Hidden = true)%Access=private
        nx
        % the elements below should be for reduced-form vars only
        % hence, the var shall inherit from the structural var
        construction_data
    end
    properties(Access=protected,Hidden = true)
        param_template
        param_to_mat_links
        all_param_names_vec
    end
    methods
        function obj=svar(varargin)
            obj=obj@rise_generic();
            if isempty(varargin)
                return
            end
            if isa(varargin{1},'svar')
                obj=varargin{1};
                varargin(1)=[];
                obj=set(obj,varargin{:});
                return
            end
            out=vartools.preliminaries(varargin{:});
            obj=rise_generic.reset(obj,out.endogenous,out.exogenous,...
                out.observables,out.markov_chains);
            quickfill={'param_template','constant','nlags','nx',...
                'construction_data'};
            for ifill=1:numel(quickfill)
                obj.(quickfill{ifill})=out.(quickfill{ifill});
            end
            % reset block exogeneity, which was already set in
            % preliminaries
            %---------------------------------------------------
            obj.endogenous.is_block_exogenous=...
                ismember(obj.endogenous.name,out.block_exogenous);
            % now build the options
            %----------------------
            obj=set(obj,'initialize');
            further_inputs=out.remains;
            obj=set(obj,further_inputs{:});
            % initialize the solution field
            %------------------------------
            obj.solution=struct();
            % link the parameters to the structural matrices and initialize
            %--------------------------------------------------------------
            [obj.param_to_mat_links,obj.parameter_values,obj.all_param_names_vec]=...
                vartools.parameters_to_matrices(...
                obj.param_template,obj.parameters.name,...
                obj.markov_chains.regimes_number);
        end
        varargout=solve(varargin)
        varargout=set_solution_to_companion(varargin)
    end
    methods(Sealed)
        function varargout=estimate(obj,varargin)
            nout=nargout;
            varargout=cell(1,nout);
            if isempty(obj)
                [varargout{1:nout}]=estimate@rise_generic(obj);
            else
                [varargout{1:nout}]=estimate@rise_generic(obj);
            end
        end
        function obj=set(obj,varargin)
            nn=length(varargin);
            % override the rise_generic set method
            if isempty(obj)|| nn==1
                obj=set@rise_generic(obj,varargin{:});
                % combine with specific elements if necessary
            else
                very_special={'vp_prior_type','estim_linear_restrictions','data'};
                redo_priors=false;
                for iobj=1:numel(obj)
                    vargs=varargin;
                    constr_data=obj(iobj).construction_data;
                    construction_data_is_modified=false;
                    options_=obj(iobj).options;
                    fr=fieldnames(constr_data);
                    fo=fieldnames(options_);
                    processed=true(1,nn);
                    for iarg=1:2:nn
                        ff=vargs{iarg};
                        val=vargs{iarg+1};
                        if any(strcmp(ff,fr))
                            if ~isequal(constr_data.(ff),val)
                                construction_data_is_modified=true;
                                constr_data.(ff)=val;
                            end
                        elseif any(strcmp(ff,fo))
                            options_.(ff)=val;
                            if any(strcmp(ff,very_special))
                                redo_priors=true;
                            end
                        else
                            processed([iarg,iarg+1])=false;
                            % don't do anything, we will use the generic method
                            % to process this later on.
                        end
                    end
                    if construction_data_is_modified
                        obj(iobj)=svar(constr_data);
                    end
                    obj.routines.likelihood=@vartools.var_likelihood;
                    obj(iobj).options=options_;
                    vargs(processed)=[];
                    if ~isempty(vargs)
                        obj(iobj)=set@rise_generic(obj(iobj),vargs{:});
                    end
                    if redo_priors && ~isempty(obj(iobj).options.data)
                        if isa(obj(iobj).options.data,'ts') && ...
                                obj(iobj).options.data.NumberOfObservations==0
                            return
                        end
                        % get the names of the estimated parameters from
                        % the information in obj.estim_param_template: the
                        % transition probabilities are automatically
                        % estimated.
                        %--------------------------------------------------
                        estim_names=create_estimated_parameters_list(obj);
                        
                        obj.estimation_restrictions=parameters_links(obj,estim_names);
                        
                        % load the data
                        %--------------
                        [obj,issue,retcode]=load_data(obj);
                        if retcode
                            error(issue)
                        end
                        obj=msvar_priors(obj,estim_names);
                    end
                end
            end
        end
    end
    methods(Static,Access=private)
        varargout=decompose_parameter(varargin)
        varargout=reformat_restriction(varargin)
    end
    methods(Access=protected)
        varargout=create_estimated_parameters_list(varargin)
    end
    methods(Hidden=true)
        varargout=simulation_engine(varargin)
        varargout=conclude_estimation(varargin)
        varargout=load_solution(varargin)
    end
    methods(Access=protected,Hidden=true)
    end
    methods(Static)
        function r=template()
            r=struct('model_class','svar',...
                'constant',true,'nlags',4,...
                'endogenous',{{}},...
                'block_exogenous',{{}},...
                'observables',{{}});
            r.markov_chains=struct('name',{},'states_expected_duration',{},'controled_parameters',{});%,'transition_matrix',{}
        end
    end
end
