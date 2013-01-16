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
     if (strcmpi(fi,'MaxIter')||strcmpi(fi,'MaxIter'))&&~isempty(options.(fi))
         bbo_options.MaxIter=options.(fi);
     elseif (strcmpi(fi,'MaxFunEvals'))&&~isempty(options.(fi)) 
         bbo_options.MaxFunEvals=options.(fi);
     elseif (strcmpi(fi,'MaxTime'))&&~isempty(options.(fi)) 
         bbo_options.MaxTime=options.(fi);
     elseif (strcmpi(fi,'MaxNodes') )&&~isempty(options.(fi))
         bbo_options.MaxNodes=options.(fi);
     end
 end
bbo_options.restrictions=PROBLEM.nonlcon;

[xx,ff,obj]=bbo_3(Objective,x0,[],lb,ub,bbo_options);

x=xx(:,1);
f=ff(1);

if obj.iterations>=obj.MaxIter || ...
        obj.funcCount>=obj.MaxFunEvals || ...
        etime(obj.finish_time,obj.start_time)>=obj.MaxTime
    exitflag=0; % disp('Too many function evaluations or iterations.')
else
    exitflag=-1; % disp('Stopped by output/plot function.')
end
        
                    
