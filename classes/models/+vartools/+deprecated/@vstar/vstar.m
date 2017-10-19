classdef vstar < generic
    % VSTAR -- vstar class
    %
    % Methods
    %
    % - [draw_parameter](vstar/draw_parameter)
    % - [estimate](vstar/estimate)
    % - [evaluate_general_restrictions](vstar/evaluate_general_restrictions)
    % - [forecast](vstar/forecast)
    % - [get](vstar/get)
    % - [growth_database](vstar/growth_database)
    % - [hessian](vstar/hessian)
    % - [initial_conditions](vstar/initial_conditions)
    % - [irf](vstar/irf)
    % - [isnan](vstar/isnan)
    % - [load_parameters](vstar/load_parameters)
    % - [log_posterior_kernel](vstar/log_posterior_kernel)
    % - [log_prior_density](vstar/log_prior_density)
    % - [mode_curvature ](vstar/mode_curvature )
    % - [plot](vstar/plot)
    % - [plot_posteriors](vstar/plot_posteriors)
    % - [plot_priors](vstar/plot_priors)
    % - [plot_priors_and_posteriors](vstar/plot_priors_and_posteriors)
    % - [posterior_sample](vstar/posterior_sample)
    % - [print_estimation_results](vstar/print_estimation_results)
    % - [pull_objective](vstar/pull_objective)
    % - [set](vstar/set)
    % - [solve](vstar/solve)
    %
    % Properties
    %
    % - [nlags](vstar/nlags)
    % - [constant](vstar/constant)
    % - [thresholds](vstar/thresholds)
    % - [deterministic](vstar/deterministic)
    % - [residuals](vstar/residuals)
    % - [transitions](vstar/transitions)
    % - [legend](vstar/legend)
    % - [endogenous](vstar/endogenous)
    % - [exogenous](vstar/exogenous)
    % - [parameters](vstar/parameters)
    % - [observables](vstar/observables)
    % - [options](vstar/options)
    % - [estimation](vstar/estimation)
    % - [solution](vstar/solution)
    % - [filtering](vstar/filtering)
    
    properties
        % number of lags
        nlags=4
        % indicator for whether the model has a constant
        constant=true
        % structure holding threshold information
        thresholds
        % structure holding information on deterministic variables other
        % than the constant
        deterministic=[]
        
    end
    
    properties(Dependent)
        % residuals of the estimation
        residuals
        % historical values for the transitions
        transitions
        
    end
    
    properties(Hidden=true)%Access=protected,
        
        num_regessors
        
        reordering_index
        
        variables_locations_in_data
        
    end
    
    methods
        
        % Constructor
        %------------
        function obj=vstar(endo_names,thresholds,nlags,constant,det_vars,varargin)
            % VSTAR -- constructor for a VSTAR object
            %
            % Syntax
            % -------
            % ::
            %
            %   obj=VSTAR(endo_names,thresholds)
            %
            %   obj=VSTAR(endo_names,thresholds,nlags)
            %
            %   obj=VSTAR(endo_names,thresholds,nlags,constant)
            %
            %   obj=VSTAR(endo_names,thresholds,nlags,constant,det_vars)
            %
            %   obj=VSTAR(endo_names,thresholds,nlags,constant,det_vars,varargin)
            %
            % Inputs
            % -------
            %
            % - **endo_names** [char|cellstr]: list of endogenous variables, with
            % possibly their descriptions in " "
            %
            % - **thresholds** [vstar_threshold]: object holding information on the
            % threshold variables and functions (see VSTAR_THRESHOLD for help)
            %
            % - **nlags** [integer|{4}]: number of lags
            %
            % - **constant** [false|{true}]: allow a constant or not
            %
            % - **det_vars** [char|cellstr|{}]: List of deterministic variables other
            % than the constant term
            %
            % - **varargin** [optional]: valid pairwise arguments for a generic RISE
            % object.
            %
            % Outputs
            % --------
            %
            % - **obj** [VSTAR]: model object
            %
            % More About
            % ------------
            %
            % Examples
            % ---------
            %
            % See also:
            
            vargs={};
            
            n=nargin;
            
            if n
                
                [endo_names,xonames,obsnames]=make_plain_names(endo_names,thresholds,det_vars);
                
                vargs={endo_names,xonames,obsnames};
                
            end
            
            obj=obj@generic(vargs{:});
            
            if n
                
                check_nargin()
                
                % deterministic variables
                %------------------------
                obj.deterministic=make_names(det_vars);
                
                if any(ismember(obj.deterministic.name,obj.endogenous.name))
                    
                    error('deterministic variables cannot be endogenous')
                    
                end
                
                obj=vstar.process_thresholds(obj,thresholds);
                
                n=length(varargin);
                
                if rem(n,2)
                    
                    error('varargin arguments must come in pairs')
                    
                end
                
                set_var()
                
                obj=set(obj,varargin{:});
                
                obj=add_to_routines(obj,'likelihood',@likelihood_smooth_transition_var);
                
            end
            
            function set_var()
                
                if ~isempty(nlags)
                    
                    if ~(isscalar(nlags) && ...
                            isnumeric(nlags) && ...
                            (floor(nlags)==ceil(nlags)) && ...
                            nlags>0)
                        
                        error('nlags should be a positive integer')
                        
                    end
                    
                    obj.nlags=nlags;
                    
                end
                
                if ~isempty(constant)
                    
                    if ~ismember(constant,[true,false])
                        
                        error('constant should be true or false')
                        
                    end
                    
                    obj.constant=constant;
                    
                end
                
            end
            
            function check_nargin()
                
                if n < 5
                    
                    det_vars={};
                    
                    if n < 4
                        
                        constant=[];
                        
                        if n < 3
                            
                            nlags=[];
                            
                            if n < 2
                                
                                error('at least two input arguments are needed')
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        varargout=solve(varargin)
        
        varargout=estimate(varargin) % overloaded
        
        function r=get.residuals(obj)
            
            [resid,~,retcode]=low_level_residuals(obj);
            
            if retcode
                
                error(decipher(retcode))
                
            end
            
            exo_names=regexp(obj.exogenous.name,'\<EPS_\w+','match');
            
            exo_names=[exo_names{:}];
            
            r=struct();
            
            for irow=1:size(resid,1)
                
                if irow==1
                    
                    prototype=ts(obj.options.estim_start_date,resid(irow,:).');
                    
                end
                
                r.(exo_names{irow})=prototype.reset_data(resid(irow,:).');
                
            end
            
        end
        
        function t=get.transitions(obj)
            
            [G,retcode]=low_level_transitions(obj);
            
            if retcode
                
                error(decipher(retcode))
                
            end
            
            t=struct();
            
            for ithresh=1:size(G,1)
                
                thresh_name=[obj.thresholds(ithresh).name,...
                    strrep(obj.thresholds(ithresh).lag,'-','_')];
                
                if ithresh==1
                    
                    prototype=ts(obj.options.estim_start_date,G(ithresh,:).');
                    
                end
                
                t.(thresh_name)=prototype.reset_data(G(ithresh,:).');
                
            end
            
        end
        
        function obj=set(obj,varargin)
            
            [obj,unprocessed]=vstar.process_options(obj,varargin{:}) ;
            
            obj=set@generic(obj,unprocessed{:});
            
        end
        
        varargout=plot(varargin)
        
    end
    
    methods(Static,Access=private)
        
        varargout=process_thresholds(varargin)
        
        varargout=set_baseline_parameters(varargin)
        
        varargout=process_options(varargin)
        
    end
    
    methods(Access=private)
        
        varargout=low_level_transitions(varargin)
        
        varargout=data_location_dispatch(varargin)
        
    end
    
    methods(Hidden = true)
        
        varargout=low_level_residuals(varargin)
        
        varargout=ball_park(varargin)
        
        varargout=load_data(varargin)
        
        function [log_post,log_lik,log_prior,Incr,retcode,obj,x1]=conclude_estimation(obj,x1)
            
            obj.estimation_under_way=false;
            
            [log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj,x1);
            
        end
        
    end
    
end

function [endo_names,xonames,obsnames]=make_plain_names(endo_names,thresholds,det_vars)

[endog]=make_names(endo_names);

endonames=endog.name;

xonames=strcat('EPS_',endonames);

xonames=put_together(xonames,xonames);

thresh_vars={thresholds.items.transition_variable};

thresh_descrpt={thresholds.items.transition_description};

[thresh_vars,thresh_descrpt]=unicize(thresh_vars,thresh_descrpt);

thresh_vars=put_together(thresh_vars,thresh_descrpt);

xonames=[xonames,det_vars,thresh_vars];

obsnames=[endo_names,det_vars,thresh_vars];

    function [x,y]=unicize(x,y)
        
        tmp=unique(x);
        
        tmp2=tmp;
        
        for ii=1:numel(tmp)
            
            bingo=find(strcmp(tmp{ii},x),1,'first');
            
            tmp2{ii}=y{bingo};
            
        end
        
        x=tmp;
        
        y=tmp2;
        
    end

    function x=put_together(x,d)
        
        d=add_dbl_qts(d);
        
        x=[x;d];
        
        x=x(:).';
        
    end

    function x=add_dbl_qts(x)
        
        x=strcat('"',x,'"');
        
    end

end

function [endog]=make_names(endo_names)

endog=cell(2,size(endo_names,2));

iter=0;

while ~isempty(endo_names)
    
    if endo_names{1}(1)=='"'
        
        if iter==0 ||...
                ~isempty(endog{2,iter})||...
                isempty(endog{1,iter})
            error(['description coming before model ',...
                'variable or two consecutive descriptions'])
            
        end
        
        endog{2,iter}=strrep(endo_names{1},'"','');
        
    else
        
        iter=iter+1;
        
        endog{1,iter}=endo_names{1};
        
    end
    
    endo_names=endo_names(2:end);
    
end

endog=endog(:,1:iter);

for ii=1:iter
    
    if isempty(endog{2,ii})
        
        endog{2,ii}=endog{1,ii};
        
    end
    
end

[~,tags]=sort(endog(1,:));

endog=endog(:,tags);

endog=struct('name',{endog(1,:)},...
    'tex_name',{endog(2,:)},'number',numel(tags));

end
