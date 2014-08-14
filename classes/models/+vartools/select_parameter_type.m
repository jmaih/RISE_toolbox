function [names,locs]=select_parameter_type(estim_names,type)
switch type
    case 'lag_coef'
        % lag matrices a0, a1,...,ap
        %---------------------------
        names=regexp(estim_names,'(?<!w+)a\d+_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        names=[names{:}];
    case 'det_coef'
        % deterministic terms
        %--------------------
        names=regexp(estim_names,'c_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        names=[names{:}];
    case 'stdev_corr'
        % standard deviations and correlations
        %-------------------------------------
        names=regexp(estim_names,'(sig|omg)_\d+_\d+(_\w+_\d+)?(?!\w+)','match');
        names=[names{:}];
    case 'theta_coef'
        % standard deviations of various processes
        %-----------------------------------------
        names=regexp(estim_names,'theta_(a\d+|c|sig|omg)_\d+_\d+(?!\w+)','match');
        names=[names{:}];
    case 'ar_coef'
        % AR coefficients on processes
        %-----------------------------
        names=regexp(estim_names,'rho_(a\d+|c|sig|omg)_\d+_\d+(?!\w+)','match');
        names=[names{:}];
    case 'trans_probs'
        % transition probabilities
        %-------------------------
        names=regexp(estim_names,'\w+_tp_\d+_\d+','match');
        names=[names{:}];
    otherwise
        error(['unknown var parameter type "',type,'"'])
end
if nargout>1
    locs=locate_variables(names,estim_names);
end
end
