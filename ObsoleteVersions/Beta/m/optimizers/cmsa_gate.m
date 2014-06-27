function [x,f,exitflag,obj]=cmsa_gate(PROBLEM)
%  CMSA_GATE attempts to find the global minimum of a constrained function of
%  several variables.
%   CMSA_GATE attempts to solve problems of the form:
%    min F(X)  subject to:  LB <= X <= UB   (bounds)
%     X
%
%   RES = CMSA_GATE(FUN,X0,LB,UB) starts at X0 and finds a minimum X to the
%   function FUN, subject to the bound constraints LB and UB. FUN accepts
%   input X and returns a scalar function value F evaluated at X. X0 may be
%   a scalar or a vector.
%
%   [X,F,OBJ] = CMSA_GATE(FUN,X0,LB,UB,OPTIONS) optimizes the function FUN under the
%   optimization options set under the structure OPTIONS. The fields of
%   this structure could be all or any of the following:
%       - 'MaxNodes': this the number of different elements in the group
%       that will share information with each other in order to find better
%       solutions. The default is 20
%       - 'MaxIter': the maximum number of iterations. The default is 1000
%       - 'MaxTime': The time budget in seconds. The default is 3600
%       - 'MaxFunEvals': the maximum number of function evaluations. The
%       default is inf
%       - 'rand_seed': the seed number for the random draws
%    X if the optimum, F is the value of the objective at X, OBJ is an
%    object with more information
%
%   Optimization stops when one of the following happens:
%   1- the number of iterations exceeds MaxIter
%   2- the number of function counts exceeds MaxFunEvals
%   3- the time elapsed exceeds MaxTime
%   4- the user write anything in and saves the automatically generated
%   file called "ManualStopping.txt"
%
%   Examples
%     FUN can be specified using @:
%        X = cmsa_gate(@myfunc,...)
%     In this case, F = myfunc(X) returns the scalar function value F of
%     the MYFUNC function evaluated at X.
%
%     FUN can also be an anonymous function:
%        X = cmsa_gate(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
%     returns X = [0;0].
%
%     FUN=inline('sum(x.^2)'); n=100;
%     lb=-20*ones(n,1); ub=-lb; x0=lb+(ub-lb).*rand(n,1);
%     optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,...
%     'MaxFunEvals',inf);
%     x=cmsa_gate(@(x) FUN(x),x0,lb,ub,optimpot)

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   $Revision: 7 $  $Date: 2011/05/26 11:23 $


Objective=PROBLEM.objective;
x0=PROBLEM.x0;
lb=PROBLEM.lb;
ub=PROBLEM.ub;
options=PROBLEM.options;

 fields=fieldnames(options);
 cmsa_options=[];
 for ii=1:numel(fields)
     fi=fields{ii};
     if (strcmpi(fi,'MaxIter')||strcmpi(fi,'MaxIter'))&&~isempty(options.(fi))
         cmsa_options.MaxIter=options.(fi);
     elseif (strcmpi(fi,'MaxFunEvals')||strcmpi(fi,'MaxFunEvals'))&&~isempty(options.(fi)) 
         cmsa_options.MaxFunEvals=options.(fi);
     elseif (strcmpi(fi,'MaxTime')||strcmpi(fi,'MaxTime'))&&~isempty(options.(fi)) 
         cmsa_options.MaxTime=options.(fi);
     elseif (strcmpi(fi,'MaxNodes')||strcmpi(fi,'MaxNodes') )&&~isempty(options.(fi))
         cmsa_options.MaxNodes=options.(fi);
     end
 end
cmsa_options.restrictions=PROBLEM.nonlcon;

obj=cmsa(Objective,x0,[],lb,ub,cmsa_options);

x=obj.best;
f=obj.best_fval;

if obj.iterations>=obj.MaxIter || ...
        obj.funcCount>=obj.MaxFunEvals || ...
        etime(obj.finish_time,obj.start_time)>=obj.MaxTime
    exitflag=0; % disp('Too many function evaluations or iterations.')
else
    exitflag=-1; % disp('Stopped by output/plot function.')
end
        
                    
