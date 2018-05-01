function stopflag=check_convergence(obj)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

stopflag=[];

if obj.iterations>=obj.MaxIter
    
    stopflag=sprintf('iterations budget of %0.0f counts exhausted',obj.MaxIter);
    
elseif etime(clock,obj.start_time)>=obj.MaxTime
    
    stopflag=sprintf('time budget of %0.0f seconds exhausted',obj.MaxTime);
    
elseif obj.funcCount>=obj.MaxFunEvals
    
    stopflag=sprintf('function evaluations budget of %0.0f counts exhausted',obj.MaxFunEvals);
    
elseif isfield(obj,'known_optimum_reached')&& obj.known_optimum_reached
    
    stopflag='known optimum reached';
    
else
    
    ms=utils.optim.manual_stopping(1);
    
    if ms==1
        
        stopflag='manually aborted (without extreme prejudice)';
        
    elseif ms==2
        
        keyboard
        
    elseif ms==3
        
        fprintf(1,['iterations (%0.0f of %0.0f), ',...
            'elapsed time (%0.0f of %0.0f),',...
            ' function evaluations (%0.0f of %0.0f) \n'],...
            obj.iterations,obj.MaxIter,...
            etime(clock,obj.start_time),obj.MaxTime,...
            obj.funcCount,obj.MaxFunEvals);
        
        pause(5)
        
    end
    
end

%if rem(obj.iterations,100)==0 || ~isempty(stopflag)
%    save([obj.optimizer,'_results'],'obj')
%end
if ~isempty(stopflag)
    
    disp(['Optimization terminated ::',upper(stopflag)])
    
end

end
