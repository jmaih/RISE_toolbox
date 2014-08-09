classdef svar < rise_generic
    properties
    end
    properties(SetAccess=protected)
        constant
        nlags
    end
    properties(SetAccess = private, Hidden = true)%Access=private
        nx
        estim_param_template        
        % the elements below should be for reduced-form vars only
		% hence, the var shall inherit from the structural var
        construction_data
    end
    properties(Access=protected,Hidden = true)
        param_template
        param_to_mat_links
    end
    properties(Constant,Access=protected)
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
			[obj.param_to_mat_links,obj.parameter_values]=...
                vartools.parameters_to_matrices(...
                obj.param_template,obj.parameters.name,...
                obj.markov_chains.regimes_number);
        end
        varargout=solve(varargin)
    end
    methods(Sealed)
        function varargout=estimate(obj,varargin)
            nout=nargout;
            varargout=cell(1,nout);
            if isempty(obj)
                [varargout{1:nout}]=estimate@rise_generic(obj);
            else
                % use the options to apply the zero restrictions on
                % individual parameters (above and beyond the earlier block
                % exogeneity
                %----------------------------------------------------------
                obj.estim_param_template=obj.param_template;
                % get the names of the estimated parameters from the
                % information in obj.estim_param_template: the transition
                % probabilities are automatically estimated.
                %----------------------------------------------------------
                estim_names=create_estimated_parameters_list(obj);
                
                obj.estimation_restrictions=parameters_links(obj,estim_names);

                % load the data
                %--------------
                [obj,issue,retcode]=load_data(obj,varargin{:});
                if retcode
                    error(issue)
                end
                obj=msvar_priors(obj,estim_names);
                if isa(obj,'stochvol')
                    [varargout{1:nout}]=posterior_simulator(obj,...
                        'mcmc_gibbs_sampler_func',@stochvol_tools.gibbs_sampler);
                else
                    obj.routines.likelihood=@vartools.var_likelihood;
                    if nout>nargout('rise_generic/estimate')
                        error('number of output arguments exceeds the maximum')
                    end
                    [varargout{1:nout}]=estimate@rise_generic(obj);
                end
            end
        end
        function obj=set(obj,varargin)
                nn=length(varargin);
            % override the rise_generic set method
            if isempty(obj)|| nn==1
                obj=set@rise_generic(obj,varargin{:});
                % combine with specific elements if necessary
            else
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
                        else
                            processed([iarg,iarg+1])=false;
                            % don't do anything, we will use the generic method
                            % to process this later on.
                        end
                    end
                    if construction_data_is_modified
                        obj(iobj)=svar(constr_data);
                    end
                    obj(iobj).options=options_;
                    vargs(processed)=[];
                    if ~isempty(vargs)
                        obj(iobj)=set@rise_generic(obj(iobj),vargs{:});
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
        varargout=load_order_var_solution(varargin)
    end
    methods(Access=private)
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
