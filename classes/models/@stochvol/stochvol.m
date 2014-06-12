classdef  stochvol < rfvar
    properties(Hidden=true)
        endogenous_positions
        exogenous_positions
        parameters_positions
    end
    properties(SetAccess=protected)
        time_varying_parameters
        random_walk_parameters
    end   
    methods
        function obj=stochvol(varargin)
            obj=obj@rfvar(varargin{:});
            if ~isempty(varargin)
                obj.time_varying_parameters=logical(varargin{1}.time_varying_parameters);
                obj.random_walk_parameters=logical(varargin{1}.random_walk_parameters);
                if ~any(obj.time_varying_parameters)
                    error('no time varying parameters and no stochastic volatility, please use the more efficient "rfvar" object instead')
                end
                % redefine endogenous,exogenous,parameters, etc.
                %-----------------------------------------------
                obj=stochvol.redo_declarations(obj);
                % link the parameters to the structural matrices
                %-----------------------------------------------
                [obj.param_to_mat_links,obj.parameter_values]=...
                    vartools.parameters_to_matrices(...
                    obj.param_template,obj.parameters.name,...
                    obj.markov_chains.regimes_number);
            end
        end
    end
    methods(Sealed)
    end
    methods(Sealed,Hidden=true)
        varargout=simulation_engine(varargin)
    end
    methods(Static,Hidden=true)
        varargout=redo_declarations(varargin)
        varargout=locate_variables_blocks(varargin)
        varargout=form_parameter_matrices(varargin)
        varargout=format_blocks(varargin)
% % % % % % % % % % % % % %         varargout=gibbs_sampler(varargin)
    end
    methods(Static)
        function r=template()
            r=rfvar.template();
            r.model_class='stochvol';
            % lags,volatility,correlations
            r.time_varying_parameters=[true,true,false];
            r.random_walk_parameters=true(1,3);
        end
    end
end

