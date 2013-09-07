function stopflag=check_convergence(obj)
stopflag=[];
if obj.iterations>=obj.MaxIter;
    stopflag='maximum number of iterations reached';
elseif etime(clock,obj.start_time)>=obj.MaxTime
    stopflag='time budget exhausted';
elseif obj.funcCount>=obj.MaxFunEvals
    stopflag='maximum number of function evaluations reached';
elseif isfield(obj,'known_optimum_reached')&& obj.known_optimum_reached
    stopflag='known optimum reached';
else
    ms=manual_stopping(1);
    if ms==1
        stopflag='manually aborted';
    elseif ms==2
        keyboard
    end
end
%if rem(obj.iterations,100)==0 || ~isempty(stopflag)
%    save([obj.optimizer,'_results'],'obj')
%end
if ~isempty(stopflag)
    disp(['Optimization terminated ::',upper(stopflag)])
end
end
