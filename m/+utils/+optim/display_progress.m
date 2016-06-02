function display_progress(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if rem(obj.iterations,obj.verbose)==0 || obj.iterations==1
    
    optimizer_nodes_nparam=sprintf('%s(#nodes=%0.0f,#params=%0.0f)',...
        obj.optimizer,obj.MaxNodes,obj.number_of_parameters);
    
    disperse=utils.optim.dispersion(obj.xx,obj.lb,obj.ub);
    
    fprintf(1,['iter: %6.0f(%6.0f)   fmin: %8.4f    stdev: %4.4f  ',...
        'f-Count:  %8.0f(%8.0f) elapsed time:  %6.0f(%6.0f)  routine: %s\n'],...
        obj.iterations,obj.MaxIter,...
        obj.best_fval,disperse,...
        obj.funcCount,obj.MaxFunEvals,...
        etime(clock,obj.start_time),obj.MaxTime,...
        optimizer_nodes_nparam);
    
end

end