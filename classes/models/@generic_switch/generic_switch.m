classdef generic_switch < generic
    
    properties(SetAccess=protected)
        % information on markov chains, regimes and related items
        markov_chains
        
    end
        
    methods(Abstract)
        % methods that must be implemented by the subclasses
        varargout=solve(varargin)
        
        % methods that must be implemented by the subclasses
        varargout=set_solution_to_companion(varargin)
        
    end
    
    methods(Abstract, Hidden = true)
        % methods that must be implemented by the subclasses
        
        varargout=problem_reduction(varargin)
        
    end
    
    methods
        
        function obj=generic_switch(varargin)
            
            if nargin
                
                obj=generic_switch.reset(obj,varargin{:});
                
            end
            
        end
        
        varargout=historical_decomposition(varargin)
        
        varargout=is_stable_system(varargin)
        
        varargout=theoretical_autocorrelations(varargin)
        
        varargout=theoretical_autocovariances(varargin)
        
        varargout=variance_decomposition(varargin)
        
    end
    
    methods(Static,Hidden=true)
        
        function obj=reset(obj,varargin)
            
            obj=reset@generic(obj,varargin(1:3));
            
            mark_parameters()
            
            markchains=varargin{4};
            
            % the markov chains will set the parameters
            %------------------------------------------
            obj=add_markov_chains_and_parameters(obj,markchains);
            
            function mark_parameters()
                
                n=obj.parameters.number;
                
                obj.parameters.is_switching=false(1,n);
                
                obj.parameters.is_trans_prob=false(1,n);
                
                obj.parameters.governing_chain=ones(1,n);
                
                % re-tag the transition probabilities
                %------------------------------------
                for iparam=1:n
                    
                    obj.parameters.is_trans_prob(iparam)=...
                        parser.is_transition_probability(obj.parameters.name{iparam});
                    
                end
                
            end
            
        end
        
    end
    
    methods(Access=private)
        
        varargout=add_markov_chains_and_parameters(varargin)
        
    end
    
    methods(Hidden=true)
        
        varargout=prepare_transition_routine(varargin)
        
    end
    
    methods(Hidden=true,Sealed)
        
        varargout=decompose_parameter_name(varargin)
    end
    
end

