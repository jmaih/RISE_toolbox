function [x,f,exitflag,H,obj]=bee_gate(Objective,x0,lb,ub,options,varargin)
% Attempts to find the global minimum of a constrained function
%
% ::
%
%    [x,f,exitflag,H,obj] = bee_gate(Objective,x0,lb,ub,options,varargin)
%
% Args:
%    Objective (function handle): objective function to minimize
%    x0 (double): initial guess of the argmin
%    lb (double): lower bound for parameters
%    ub (double): upper bound for parameters
%    options (struct): optimize options
%
%       - **'MaxNodes'**:{20} this the number of different elements in the group
%         that will share information with each other in order to find better
%         solutions.
%       - **'MaxIter'**: {1000} the maximum number of iterations.
%       - **'MaxTime'**: {3600} The time budget in seconds.
%       - **'MaxFunEvals'**: {inf} the maximum number of function evaluations.
%       - **'rand_seed'**: the seed number for the random draws (for replicability)
%
%    varargin : extra arguments to be fed into the function to minimize
%
% Returns:
%    :
%
%    - x (double): argmin of the Objective function
%    - f (double): minimum value found by csminwel
%    - exitflag (int): **Currently(2018/07/12) turned off**
%
%       - 1: the number of iterations exceeds MaxIter
%       - 2: the number of function counts exceeds MaxFunEvals
%       - 3: the time elapsed exceeds MaxTime
%       - 4: the user write anything in and saves the automatically generated
%         file called "ManualStopping.txt"
%
%    - H (matrix): EXPERIMENTAL! Inverse of the covariance of the best values (if unimodal, can be used as an approximation of the covariance)
%    - obj (bee object): INTERNAL INFO! struct containing all information from the optimization process. Support for this class is not provided
%
% Example:
%    ::
%
%       % Simple Example
%       X = bee_gate(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0]);
%
%       % More Elaborate
%       FUN=inline('sum(x.^2)');
%       n=100;
%       lb=-20*ones(n,1);
%       ub=-lb;
%       x0=lb+(ub-lb).*rand(n,1);
%       optimpot=struct('MaxNodes',20,'MaxIter',1000,'MaxTime',60,'MaxFunEvals',inf);
%
%       x=bee_gate(@(x) FUN(x),x0,lb,ub,optimpot)
%
% References:
%     - :cite:`karaboga2014comprehensive`
%

%   Copyright 2011 Junior Maih (junior.maih@gmail.com).
%   Revision: 7
%   Date: 2011/05/26 11:23

% OtherProblemFields={'Aineq','bineq','Aeq','beq','nonlcon'};


bee_fields= properties('bee');

bee_options=struct();

for ii=1:numel(bee_fields)

    fi=bee_fields{ii};

    if isfield(options,fi)

        bee_options.(fi)=options.(fi);

    end

end

f0=[];

if iscell(x0)

    if numel(x0)~=2

        error('when x0 is a cell, it should have two elements')

    end

    f0=x0{2};

    x0=x0{1};

end

% if license('checkout','Distrib_Computing_Toolbox') && matlabpool('size') &&...
% 	(isfield(options,'UseParallel') && strcmp(options.UseParallel,'always'))
%     obj=par_bee(Objective,x0,[],lb,ub,bee_options,varargin{:});
% else
    obj=bee(Objective,x0,f0,lb,ub,bee_options,varargin{:});
% end

x=obj.best;

f=obj.best_fval;

nx=numel(x);

H=obj.COV\eye(nx);

if obj.iterations>=obj.MaxIter || ...
        obj.funcCount>=obj.MaxFunEvals || ...
        etime(obj.finish_time,obj.start_time)>=obj.MaxTime

    exitflag=0; % disp('Too many function evaluations or iterations.')

else

    exitflag=-1; % disp('Stopped by output/plot function.')

end

end
