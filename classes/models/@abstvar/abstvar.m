classdef (Abstract) abstvar < gogetter
    
    properties(SetAccess=protected)
        % you cannot set those
        endogenous
        exogenous
        parameters
        members
        linear_restrictions
        data
        constant
        homogeneity
        debug=false
        markov_chains
        date_range
        estim_param
        prior 
    end
    
    properties(SetAccess=protected,Hidden)
        
        linear_restrictions_prime_time = cell(1,0)
        
        markov_chain_info
        
    end
    
    properties(SetAccess=protected,Hidden)
        % you cannot set those and they are hidden
        % kdata
        nlags
        nx
        ng
        nvars
        Y
        X
        K % number of estimated parameters per equation
        T % number of observations after treating the lags
        sampler
        mapping
        linres
        nonlinres
    end
    
    properties(Dependent)
        is_panel
        nparams
        nregs
        is_switching
    end
    
    properties(Constant)
        
        panel_types={'pooled','meanGroup','static','dynamic','independent','unrestricted'}
        
        panel_with_constant={'pooled','meanGroup','static','dynamic','independent'}
        % - pooled: everything pooled: = dynamic + static
        % - meanGroup: average across groups
        % - static: static homogeneity: deterministic coefficients common
        % - dynamic: dynamic homogeneity: lag coefficients common,
        %   different constants
        % - unrestricted: no constant
        % - independent: individual VARs
        
    end
    
    methods(Access=protected)
        
        varargout=setOptions(varargin)
        
        varargout=collect_data(varargin)
        
        varargout=set_data_to_time_series(varargin)
        
        varargout=process_linear_restrictions(varargin)
        
    end
    
    methods(Access=protected,Hidden)
        
        varargout=prime_time(varargin)
        
        varargout=stretch_variables(varargin)
        
        varargout=print_low_level(varargin)
        
        varargout=setup_nonlinear_restrictions(varargin)
        
    end
    
    methods(Static,Access=private)
                
    end
    
    methods(Static)
        
        function prior=prior_template()
            % L1 : Overall tightness
            % L2 : Cross-variable specific variance parameter
            % L3 : Speed at which lags greater than 1 converge to zero
            % L4 : tightness on deterministic/exogenous terms
            % L5 : covariance dummies(omega)
            % L6 : co-persistence (lambda)
            % L7 : Own-persistence (mu)
            
            prior=vartools.prior_hyperparams();
            
        end
        
    end
    
    methods(Static,Hidden)
        
        varargout=embed(varargin)
        
        varargout=create_variable_names(varargin)
        
        varargout=reset_range(varargin)
        
        varargout=create_parameters_names(varargin)
        
        varargout=map_estimation(varargin)
        
        varargout=recreate_parameters(varargin)
        
        varargout=map_panel(varargin)
        
        varargout=reinflate(varargin)
                
    end
    
    methods
        
        function self=abstvar(endog,exog,nlags,const,panel,markov_chains)%,varargin
            
            self=self@gogetter();
            
            n=nargin;
            
            if n>0
                
                set_defaults()
                
                check_inputs()
                
                self.endogenous=endog(:)';
                
                self.exogenous=exog(:)';
                
                self.members=panel.members(:)';
                
                self.homogeneity=panel.homogeneity;
                
                self.markov_chains=markov_chains;
                
                self.constant=const;
                
                ng=numel(panel.members);
                
                ng=max(1,ng);
                
                nx=numel(exog)+self.constant;
                
                self.nlags=nlags;
                
                self.nx=nx;
                
                self.ng=ng;
                
                self.nvars=numel(endog);
                
            end
            
            function set_defaults()
                
                if n<6
                    
                    markov_chains=[];
                    
                    if n<5
                        
                        panel=[];
                        
                        if n<4
                            
                            const=[];
                            
                            if n<3
                                
                                nlags=[];
                                
                                if n<2
                                    
                                    exog=[];
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                if isempty(const),const=true; end
                
                if isempty(nlags),nlags=4; end
                
                if isempty(panel),panel=cell(1,0); end
                
                if isempty(exog),exog=cell(1,0); end
                
            end
            
            function check_inputs()
                
                endog=cellstring_check(endog,'endogenous',false);
                
                exog=cellstring_check(exog,'exogenous',true);
                
                badvars=intersect(exog,endog);
                
                if ~isempty(badvars)
                    
                    disp(badvars)
                    
                    error('The variables above have been declared as both endogenous and exogenous')
                    
                end
                
                if isempty(panel)
                    
                    panel=struct('members',{cell(1,0)},'homogeneity','pooled');
                    
                end
                
                if ~isfield(panel,'homogeneity')
                    
                    panel.homogeneity='pooled';
                    
                end
                    
                check_lags()
                
                check_constant(const,'constant (4th input)')
                
                check_markov_chains()
                
                check_panel()
                
                function check_panel()
                    
                    if ~isstruct(panel)||~isfield(panel,'members')
                        
                        error(['panel should be a struct with fields members ',...
                            'and optionally homogeneity'])
                        
                    end
                    
                    ff=fieldnames(panel);
                    
                    bad=~ismember(ff,{'members','homogeneity'});
                    
                    if any(bad)
                        
                        disp(ff(bad))
                        
                        error('The fields above are not valid for the panel structure')
                        
                    end
                    
                    
                    if ~(ischar(panel.homogeneity) && ...
                            ismember(panel.homogeneity,self.panel_types))
                        
                        disp(self.panel_types)
                        
                        error('homogeneity must be one of the above')
                        
                    end
                    
                    panel.members=cellstring_check(panel.members,...
                        'panel.members',true);
                    
                    if numel(panel.members) && ~ismember(panel.homogeneity,...
                            self.panel_with_constant)
                        
                        if const
                            
                            warning(['not possible to have a constant in a "',...
                                panel.homogeneity,'" panel. Changing constant to false'])
                            
                            const=false;
                            
                        end
                        
                    end
                    
                end
                
                function check_markov_chains()
                    
                    if isempty(markov_chains)
                        
                        return
                        
                    end
                    
                    expected_fields=sort({'name','states_expected_duration',...
                        'controlled_parameters'});
                    
                    ff=sort(fieldnames(markov_chains));
                    
                    if ~(all(ismember(expected_fields,ff)) && ...
                            all(ismember(ff,expected_fields)))
                        
                        disp(expected_fields)
                        
                        error('markov chains can only contain the fields above')
                        
                    end
                    
                    mcNames={markov_chains.name};
                    
                    if numel(unique(mcNames))~=numel(mcNames)
                        
                        error('several markov chains with the same name')
                        
                    end
                    
                end
                
                function x=cellstring_check(x,inpName,emptyPossible)
                    
                    if isempty(x)
                        
                        if emptyPossible
                            
                            return
                            
                        end
                        
                        error([inpName,' cannot be empty'])
                        
                    end
                    
                    x=vartools.cellstringize(x);
                        
                end
                
                function check_lags()
                    
                    if isempty(nlags)
                        
                        error('nlags (3rd input) cannot be empty')
                        
                    end
                    
                    if ~(isscalar(nlags) &&...
                            isnumeric(nlags) &&...
                            isreal(nlags) &&...
                            nlags >0 && ...
                            isfinite(nlags) &&...
                            floor(nlags)==ceil(nlags))
                        
                        error('nlags (3rd input) must be a finite positive integer')
                        
                    end
                    
                end
                
                function check_constant(const,msgIntro)
                    
                    if ~(isscalar(const) && islogical(const))
                        
                        error([msgIntro,' must be true or false'])
                        
                    end
                    
                end
                
            end
            
        end
        
        function n=get.nparams(self)
            
            n=numel(self.parameters);
            
        end
        
        function n=get.nregs(self)
            
            n=size(self.mapping.regimes,1)-1;
            
        end
        
        function ok=get.is_panel(self)
            
            ok=self.ng>1;
            
        end
        
        function ok=get.is_switching(self)
            
            ok=numel(self.markov_chains)>1;
            
        end
                
        varargout=filter(varargin)
                
        varargout=solve(varargin)
                
        varargout=autocov(varargin)
                
        varargout=estimate(varargin)
        
        varargout=print_solution(varargin)
                
        varargout=residuals(varargin)
        
        varargout=set(varargin)
        
        varargout=variance_decomposition(varargin)
        
        varargout=historical_decomposition(varargin)
        
        varargout=irf(varargin)
        
        varargout=forecast(varargin)
        
    end
    
end
