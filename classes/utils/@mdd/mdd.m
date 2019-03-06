classdef mdd<handle
    
    properties
        
        log_post_kern
        
        theta_draws
        
        lb
        
        ub
        
        maximization=false
        
    end
    
    properties(Hidden)
        
        facility=true
        
        facilitator
        
        d
        
        M
        
        theta_mode
        
        recenter
        
        moms
        
        LogPost_M
        
    end
    
    properties(Dependent)
        
        thecoef
                
    end
    
    methods
        
        function obj=mdd(theta_draws,log_post_kern,lb,ub,maximization)
            
            n=nargin;
            
            if n
                
                set_inputs()
                
            end
            
            function set_inputs()
                
                if n>1
                    
                    assert(isempty(log_post_kern)||...
                        isa(log_post_kern,'function_handle'),...
                        'log_posterior_kern should be empty or a function handle')
                    
                    obj.log_post_kern=log_post_kern;
                    
                    if n>2
                        
                        if n<3
                            
                            error('wrong number of arguments')
                            
                        end
                        
                        obj.lb=lb;
                        
                        obj.ub=ub;
                        
                        if n>4
                            
                            assert(islogical(maximization) && ...
                                isscalar(maximization),...
                                'maximization should be a logical scalar')
                            
                            obj.maximization=maximization;
                            
                        end
                        
                    end
                    
                end
                
                [theta_draws,fd]=mcmc.process_draws(theta_draws);

				theta_draws=theta_draws(:,:);
				
				fd=fd(:,:);
                
                c=obj.thecoef;
                
                % multiply accordingly
                obj.LogPost_M=c*fd;
                
                [~,best_loc]=max(obj.LogPost_M);
                
                % the exponential explodes very quickly. We apply the following fix
                %-------------------------------------------------------------------
                obj.facilitator=obj.facility*obj.LogPost_M(best_loc);
                
                obj.theta_mode=theta_draws(best_loc(1));
                
                obj.theta_draws=theta_draws;  clear theta_draws
                
                [obj.d,obj.M] = size(obj.theta_draws);
                
                if isempty(obj.lb)
                    
                    obj.lb=min(obj.theta_draws,[],2);
                    
                end
                
                if ~isequal(size(obj.lb),[obj.d,1])
                    
                    error('wrong format of lower bound')
                    
                end
                
                if isempty(obj.ub)
                    
                    obj.ub=max(obj.theta_draws,[],2);
                    
                end
                
                if ~isequal(size(obj.ub),[obj.d,1])
                    
                    error('wrong format of upper bound')
                    
                end
                
                
                % recentering function
                %----------------------
                obj.recenter=@(x)utils.optim.recenter(x,obj.lb,obj.ub,3);
                
            end
            
        end
        
        function c=get.thecoef(obj)
            
            c=2*obj.maximization-1;
            
        end
               
        varargout=mhm(obj,varargin)
        
        varargout=mueller(obj,varargin)
        
        varargout=bridge(obj,varargin)
        
        varargout=is(obj,varargin)
        
        varargout=ris(obj,varargin)
        
        varargout=cj(obj,varargin)
        
        varargout=swz(obj,varargin)
        
        varargout=laplace(obj,varargin)
        
        varargout=laplace_mcmc(obj,varargin)
              
    end
    
    methods(Access=private)
        
        varargout=moments(varargin)
        
        varargout=iid_draws(varargin)
        
        varargout=old_draws_in_weighting_function(varargin)
        
    end
    
    methods(Static)
        
        varargout=normal_weighting(varargin)
        
        varargout=global_options(varargin)
        
    end
    
end