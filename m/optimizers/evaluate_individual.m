function indiv=evaluate_individual(bird,objective,lb,ub,nonlcon,varargin)
if nargin<1
    indiv=struct('x',{},'f',{},'viol',{},'violstrength',{});
else
    % correct the bounds
    bird(bird<lb)=lb(bird<lb);
    bird(bird>ub)=ub(bird>ub);
    % evaluate
    fval=objective(bird,varargin{:});
    % evaluate the constraints if any
    % I think a persistent variable that recognizes the input could
    % work well but that is not the concern of the optimizer
    if isempty(nonlcon)
        vv=0;
    else
        vv=nonlcon(bird);
    end
    indiv=struct('x',bird,...
        'f',fval,...
        'viol',vv,...
        'violstrength',sum(vv(vv>0)));
end
end
