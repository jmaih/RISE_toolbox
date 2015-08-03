function y=auxiliary_endo_sstate_evaluation(obj,y,xss,param,def)
% auxiliary_endo_sstate_evaluation -- computes the steady state for the
% auxiliary variables created by RISE
%
% Syntax
% -------
% ::
%
%   y=auxiliary_endo_sstate_evaluation(obj,y,xss,param,def)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object
%
% - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
% - **xss** [vector]: exo_nbr x 1 vector of shocks steady state values
%
% - **param** [vector]: vector of parameter values
%
% - **def** [vector]: vector of definitions
%
% Outputs
% --------
%
% - **y** [vector]: endo_nbr x 1 vector of updated steady state
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

aux_ssfunc=obj.routines.steady_state_auxiliary_eqtns;

ss=y;
if ~isempty(aux_ssfunc) && (isa(aux_ssfunc,'function_handle')||...
        (isstruct(aux_ssfunc)&&~isempty(aux_ssfunc.code)))
    % y, x, ss, param, def, s0, s1
    y=utils.code.evaluate_functions(aux_ssfunc,y,xss,ss,param,def,[],[]);
end

if obj.is_optimal_policy_model
    siz=obj.routines.planner_static_mult_support{1};
    pos=obj.routines.planner_static_mult_support{2};
    
    vals=zeros(size(pos));
    vals(pos)=utils.code.evaluate_functions(obj.routines.planner_static_mult,...
        y,xss,ss,param,def,[],[]); % y, x, ss, param, def, s0, s1
    vals=reshape(vals,siz);
    
    wx=vals(end,:);
    vals(end,:)=[];
    
    multvals=vals.'\wx.';
    
    y(obj.endogenous.is_lagrange_multiplier)=multvals;
end

end