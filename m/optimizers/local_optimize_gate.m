function [x,f,exitflag]=local_optimize_gate(PROBLEM)

Objective=PROBLEM.objective;
x0=PROBLEM.x0;
lb=PROBLEM.lb;
ub=PROBLEM.ub;
options=PROBLEM.options;

 fields=fieldnames(options);
 new_options=[];
 for ii=1:numel(fields)
     fi=fields{ii};
     if (strcmpi(fi,'max_iter')||strcmpi(fi,'MaxIter'))&&~isempty(options.(fi))
         new_options.max_iter=options.(fi);
     elseif (strcmpi(fi,'MaxFunEvals')||strcmpi(fi,'max_fcount'))&&~isempty(options.(fi)) 
         new_options.max_fcount=options.(fi);
     elseif (strcmpi(fi,'max_time')||strcmpi(fi,'MaxTime'))&&~isempty(options.(fi)) 
         new_options.max_time=options.(fi);
     elseif (strcmpi(fi,'MaxNodes')||strcmpi(fi,'colony_size') )&&~isempty(options.(fi))
         new_options.colony_size=options.(fi);
     end
 end
new_options.restrictions=PROBLEM.nonlcon;

[x,f]=local_optimize(Objective,x0,[],lb,ub,new_options);

%if obj.iter>=obj.max_iter || ...
%        obj.fcount>=obj.max_fcount || ...
%        etime(obj.finish_time,obj.start_time)>=obj.max_time
%    exitflag=0; % disp('Too many function evaluations or iterations.')
%else
    exitflag=-1; % disp('Stopped by output/plot function.')
%end
        
                    
