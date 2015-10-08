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
    % print_estimation_results -   H1 line
    % prior_plots -   H1 line
    % refresh -  refresh the options of an old object with a newer version of
    % report - assigns the elements of interest to a rise_report.report object
    % set -  sets options for RISE models
    % set_solution_to_companion -   H1 line
    % simulate -  simulates a RISE model
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
            r=varargin{1};
            % make sure you do not mess with the prior hyperparameters
            %---------------------------------------------------------
            svar.check_template(r,class(obj))
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
        varargout=filter(varargin)
        varargout=solve(varargin)
        varargout=set_solution_to_companion(varargin)
    end
    methods(Sealed)
        varargout=estimate(varargin)
        function obj=set(obj,varargin)
            nn=length(varargin);
            % override the rise_generic set method
            if isempty(obj)|| nn==1
                obj=set@rise_generic(obj,varargin{:});
                % combine with specific elements if necessary
            else
                very_special={'estim_linear_restrictions','data'};
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
                        
                        % load the data
                        %--------------
                        [obj,issue,retcode]=load_data(obj);
                        if retcode
                            error(issue)
                        end
						% msvar_priors calls setup_priors, which does parameter_links
                        obj=msvar_priors(obj);
                    end
                end
            end
        end
    end
    methods(Static,Access=private)
        function out=check_template(r,model_class)
            initiate=nargin==0;
            if initiate
                r=struct('endogenous',{{}});
            end
            num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);
            mycellstr=@(x)(numel(cellstr(x))==1 && ismember(cellstr(x),{'all','none'}))||...
                all(ismember(x,r.endogenous));
            prior_types={'minnesota','none','normal_wishart',...
                'indep_normal_wishart','jeffrey','diffuse'};
            mydefaults={
                'vp_mnst_overall_tightness',3,@(x)num_fin(x) && x>0,...
                'vp_mnst_overall_tightness must be postive and finite'
                
                'vp_mnst_relative_tightness_lags',1,@(x)num_fin(x) && x>0,...
                'vp_mnst_relative_tightness_lags must be postive and finite'
                
                'vp_mnst_relative_tightness_constant',0.1,@(x)num_fin(x) && x>0,...
                'vp_mnst_relative_tightness_lags must be postive and finite'
                
                'vp_mnst_tightness_on_lag_decay',0.5,@(x)num_fin(x) && x>0,...
                'vp_mnst_relative_tightness_lags must be postive and finite'
                
                'vp_mnst_unit_root_vars','all',@(x)mycellstr(x),...
                'vp_mnst_unit_root_vars must be "all", "none", a char or a cellstr'
                
                'vp_mnst_stationary_var_mean',0.5,@(x)num_fin(x) && x>=0 && x<1,...
                'vp_mnst_stationary_var_mean must be in [0,1)'
                
                'vp_natconj_normwish_variance',10,@(x)num_fin(x) && x>0,...
                'vp_natconj_normwish_variance must be postive and finite'
                
                % use the covariance formed by the ar1 processes when forming the posterior mode( else use the covariance of the OLS)
                'vp_gls_ar1_processes',true,@(x)islogical(x),...
                'vp_gls_ar1_processes must be a logical' 
                
                'vp_prior_type','minnesota',@(x)ismember(x,prior_types),...
                ['vp_prior_type must be one of ',cell2mat(strcat(prior_types,'|'))]
                };
            if initiate
                flag=cell2struct(mydefaults(:,2),mydefaults(:,1),1);
            else
                if ~strcmp(r.model_class,model_class)
                    error(['model object''s class (',model_class,') ',...
                        'different from template class (',r.model_class,')'])
                end
                flag=r;
                flag.priors_hyperparams=parse_arguments(mydefaults,r.priors_hyperparams);
            end
            if nargout
                out=flag;
            end
        end
        varargout=decompose_parameter(varargin)
        varargout=reformat_restriction(varargin)
    end
    methods(Access=protected)
        varargout=create_estimated_parameters_list(varargin)
    end
    methods(Hidden=true)
        varargout=simulation_engine(varargin)
        varargout=conclude_estimation(varargin) % abstract method
        varargout=load_solution(varargin)
        varargout=problem_reduction(varargin) % abstract method
    end
    methods(Access=protected,Hidden=true)
    end
    methods(Static)
        function r=template()
            % TEMPLATE - creates a template for SVAR models
            %
            % Syntax
            % -------
            % ::
            %
            %   r=svar.template()
            %
            % Inputs
            % -------
            %
            % none
            %
            % Outputs
            % --------
            %
            % - **r** [struct]: structure that can be modified by the user
            % before creating the SVAR system.
            %   - **model_class** [char|{'svar'}]: class of the model
            %   - **constant** [false|{true}]: decides whether or not the
            %   model includes a constant term.
            %   - **nlags** [integer|{4}]: number of lags
            %   - **endogenous** [cellstr|{}]: place holder of endogenous
            %   variables. The names are entered as a cellstr but they can
            %   also contain the description of the variables, in which
            %   case each model variable name is followed by its
            %   description in double quotes. e.g. {'C','"Consumption"'}.
            %   The names can also be entered as two-column cell array, in
            %   which case the first column is the model variable names and
            %   the second column is their description.
            %   - **block_exogenous** [cellstr|{}]: list of variables that
            %   are block exogenous
            %   - **observables** [cellstr|{}]: list of deterministic
            %   variables other than the constant term
            %   - **priors_hyperparams** [struct]:
            %       - **vp_mnst_overall_tightness** [numeric|{3}]:
            %       - **vp_mnst_relative_tightness_lags** [numeric|{1}]:
            %       - **vp_mnst_relative_tightness_constant** [|{0.1}]:
            %       - **vp_mnst_tightness_on_lag_decay** [|{0.5}]:
            %       - **vp_mnst_unit_root_vars** [char|cellstr|'none'|{'all'}]:
            %       - **vp_mnst_stationary_var_mean** [numeric|{0.5}]:
            %       - **vp_natconj_normwish_variance** [numeric|{10}]:
            %       - **vp_gls_ar1_processes** [false|{true}]:
            %       - **vp_prior_type** ['none' | 'normal_wishart' |
            %       'indep_normal_wishart' | 'jeffrey' | 'diffuse' |
            %       {'minnesota'}] : type of prior
            %   - **markov_chains** [struct]:
            %       - **name** [char|{}]:
            %       - **states_expected_duration** [row vector|{}]:
            %       - **controled_parameters** [cellstr|{}]:
            %
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
            priors=svar.check_template();
            r=struct('model_class','svar',...
                'constant',true,'nlags',4,...
                'endogenous',{{}},...
                'block_exogenous',{{}},...
                'observables',{{}},...
                'priors_hyperparams',priors);
            r.markov_chains=struct('name',{},...
                'states_expected_duration',{},...
                'controled_parameters',{});
        end
    end
end
