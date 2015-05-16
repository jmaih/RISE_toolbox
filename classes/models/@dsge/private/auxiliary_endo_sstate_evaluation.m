function y=auxiliary_endo_sstate_evaluation(obj,y,xss,param,def)

aux_ssfunc=obj.routines.steady_state_auxiliary_eqtns;

if ~isempty(aux_ssfunc) && (isa(aux_ssfunc,'function_handle')||...
        (isstruct(aux_ssfunc)&&~isempty(aux_ssfunc.code)))
    % y, x, ss, param, sparam, def, s0, s1
    y=utils.code.evaluate_functions(aux_ssfunc,y,xss,[],param,[],[],[],[]);
end

if obj.is_optimal_policy_model
    siz=obj.routines.planner_static_mult_support{1};
    pos=obj.routines.planner_static_mult_support{2};
    
    vals=zeros(size(pos));
    vals(pos)=utils.code.evaluate_functions(obj.routines.planner_static_mult,...
        y,xss,[],param,[],def,[],[]); % y, x, ss, param, sparam, def, s0, s1
    vals=reshape(vals,siz);
    
    wx=vals(end,:);
    vals(end,:)=[];
    
    multvals=vals'\wx';
    
    y(obj.endogenous.is_lagrange_multiplier)=multvals;
end

end