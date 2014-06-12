classdef rfvar < svar
    properties
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
    end
    methods(Sealed)
    end
    methods(Access=private)
        varargout=translate_restrictions(varargin)
        varargout=set_structural_shocks(varargin)
    end
    methods(Static,Access=private)
        varargout=sort_Q(varargin)
        varargout=var_rotation(varargin)
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

