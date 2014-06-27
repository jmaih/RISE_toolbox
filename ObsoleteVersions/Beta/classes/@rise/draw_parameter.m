function param=draw_parameter(obj,type)
if isempty(obj)
    param=struct();
    return
end
if nargin<2
    type='uniform';
end
temporary=true;
lb=vertcat(obj.estimated_parameters.lb);
ub=vertcat(obj.estimated_parameters.ub);
if temporary
    disp([mfilename,':: temporary thing should be removed'])
    if ismember(type,{'prior','prior_density'})
        type='uniform';
    end
end
switch type
    case {'prior','prior_density'}
        error([mfilename,':: prior draws not implemented yet'])
    case {'posterior','posterior_density'}
        error([mfilename,':: posterior draws not implemented yet'])
    case {'mode','mode_density'}
        C=obj.vcov;
        if temporary
            C=diag(sqrt(abs(diag(C))));
        else
            C=chol(C,'lower');
        end
        m=vertcat(obj.estimated_parameters.mode);
        notdone=true;
        while notdone
            param=m+C*randn(numel(m),1);
            if all(param>=lb) && all(param<=ub)
                notdone=false;
            end
        end
    case {'uniform'}
        param=lb+(ub-lb).*rand(numel(lb),1);
    otherwise
end
end