function stopflag=check_convergence(obj)
stopflag=[];
if obj.iterations>=obj.MaxIter;
    stopflag='maximum number of iterations reached';
end
if etime(clock,obj.start_time)>=obj.MaxTime
    stopflag='time budget exhausted';
end
if obj.funcCount>=obj.MaxFunEvals
    stopflag='maximum number of function evaluations reached';
end
if manual_stopping(1)
    stopflag='manually aborted';
end
if isfield(obj,'known_optimum_reached')&& obj.known_optimum_reached
    stopflag='known optimum reached';
end
%if rem(obj.iterations,100)==0 || ~isempty(stopflag)
%    save([obj.optimizer,'_results'],'obj')
%end
if ~isempty(stopflag)
    disp(['Optimization terminated ::',upper(stopflag)])
end
end
