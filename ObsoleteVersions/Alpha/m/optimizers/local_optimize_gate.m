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
     if (strcmpi(fi,'MaxIter'))&&~isempty(options.(fi))
         new_options.MaxIter=options.(fi);
     elseif (strcmpi(fi,'MaxFunEvals'))&&~isempty(options.(fi)) 
         new_options.MaxFunEvals=options.(fi);
     elseif (strcmpi(fi,'MaxTime'))&&~isempty(options.(fi)) 
         new_options.MaxTime=options.(fi);
     elseif (strcmpi(fi,'MaxNodes') )&&~isempty(options.(fi))
         new_options.MaxNodes=options.(fi);
     end
 end
new_options.restrictions=PROBLEM.nonlcon;

[x,f]=local_optimize(Objective,x0,[],lb,ub,new_options);

%if obj.iterations>=obj.MaxIter || ...
%        obj.funcCount>=obj.MaxFunEvals || ...
%        etime(obj.finish_time,obj.start_time)>=obj.MaxTime
%    exitflag=0; % disp('Too many function evaluations or iterations.')
%else
    exitflag=-1; % disp('Stopped by output/plot function.')
%end
        
                    
