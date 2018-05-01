function display_progress(obj)
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

if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
    
    optimizer_nodes_nparam=sprintf('%s(#nodes=%0.0f,#params=%0.0f)',...
        obj.optimizer,obj.MaxNodes,obj.number_of_parameters);
    
    fprintf(1,['iter:  %6.0f(%6.0f) ',...
        'f-Count:  %4.0f(%4.0f) elapsed time:%6.0f(%6.0f)  routine: %s\n'],...
        obj.iterations,obj.MaxIter,...
         obj.funcCount,obj.MaxFunEvals,...
        etime(clock,obj.start_time),obj.MaxTime,...
        optimizer_nodes_nparam);
    
    if isfield(obj,'accept_ratio')
        
        fprintf(1,'\t\tacceptance rates: %s\n',num2str(100*obj.accept_ratio));
        
    elseif isfield(obj,'xx')
        
        disperse=utils.optim.dispersion(obj.xx,obj.lb,obj.ub);
        
        fprintf(1,'\t\tstdev: %4.4f\n',disperse);
        
    end
    
    fprintf(1,'\t\tfmin: %s\n\n',num2str(obj.best_fval));
    
end

end