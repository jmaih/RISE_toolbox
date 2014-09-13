classdef rfvar < svar
    properties(Hidden = true)
        constant_var_data
    end
    properties(SetAccess=protected)
        identification
% % %         structural_shocks
        nonlinear_restrictions
    end
    properties(Access=protected, Hidden = true)
        % the elements below should be for reduced-form vars only
		% hence, the var shall inherit from the structural var
        permute_irf_restrictions=[2,1,3]
        permute_lag_structure_restrictions=[1,2,3]
    end
   
    methods
        function obj=rfvar(varargin)
            obj=obj@svar(varargin{:});
        end
        varargout=structural_form(varargin)
        varargout=check_identification(varargin)
        varargout=solve(varargin)
    end
    methods(Access=private)
        varargout=translate_restrictions(varargin)
        varargout=set_structural_shocks(varargin)
    end
    methods(Static,Access=private)
        varargout=sort_Q(varargin)
        varargout=var_rotation(varargin)
    end
    methods(Access=protected,Hidden=true)
        varargout=update_estimated_parameter_names(varargin)
        varargout=find_posterior_mode(varargin)
        varargout=update_posterior_simulation_initial_conditions(varargin)
        varargout=initialize_posterior_simulation(varargin)
    end
    methods(Static)
        function r=template()
            r=svar.template();
            r.model_class='rfvar';
%             r.irf_sign_restrictions=[];
%             r.irf_zero_restrictions=[];
        end
    end
end

