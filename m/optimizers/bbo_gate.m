function [x,f,exitflag,obj]=bbo_gate(PROBLEM)

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   $Revision: 7 $  $Date: 2011/05/26 11:23 $

% OtherProblemFields={'Aineq','bineq','Aeq','beq','nonlcon'};
Objective=PROBLEM.objective;
x0=PROBLEM.x0;
lb=PROBLEM.lb;
ub=PROBLEM.ub;
options=PROBLEM.options;

 fields=fieldnames(options);
 bbo_options=[];
 for ii=1:numel(fields)
     fi=fields{ii};
     if (strcmpi(fi,'max_iter')||strcmpi(fi,'MaxIter'))&&~isempty(options.(fi))
         bbo_options.max_iter=options.(fi);
     elseif (strcmpi(fi,'MaxFunEvals')||strcmpi(fi,'max_fcount'))&&~isempty(options.(fi)) 
         bbo_options.max_fcount=options.(fi);
     elseif (strcmpi(fi,'max_time')||strcmpi(fi,'MaxTime'))&&~isempty(options.(fi)) 
         bbo_options.max_time=options.(fi);
     elseif (strcmpi(fi,'MaxNodes')||strcmpi(fi,'colony_size') )&&~isempty(options.(fi))
         bbo_options.colony_size=options.(fi);
     end
 end
bbo_options.restrictions=PROBLEM.nonlcon;

[xx,ff,obj]=bbo_3(Objective,x0,[],lb,ub,bbo_options);

x=xx(:,1);
f=ff(1);

if obj.iter>=obj.max_iter || ...
        obj.fcount>=obj.max_fcount || ...
        etime(obj.finish_time,obj.start_time)>=obj.max_time
    exitflag=0; % disp('Too many function evaluations or iterations.')
else
    exitflag=-1; % disp('Stopped by output/plot function.')
end
        
                    
