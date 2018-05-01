function [x0,f0,exitflag,output]=blockwise_optimization(func,x0,lb,ub,options,varargin)
% blockwise_optimization optimization by blocks of parameters rather than the whole vector
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


if nargin<4
    options=[];
end
if isstruct(func)
    x0=func.x0;
    lb=func.lb;
    ub=func.ub;
    options=func.options;
    if isfield(func,'solver')
        options.optimizer=func.solver;
    end
    func=func.objective;
end
if ischar(func)
    func=str2func(func);
end

default_options=struct('MaxNodes',20,'MaxIter',1000,...
    'MaxTime',3600,'MaxFunEvals',inf,...
    'penalty',1e+8,'verbose',10,'blocks',[],'optimizer',@bee_gate,...
    'Display','iter','nonlcon',[]);
ff=fieldnames(default_options);
for ifield=1:numel(ff)
    v=ff{ifield};
    if isfield(options,v)
        default_options.(v)=options.(v);
    end
end

blocks=default_options.blocks;
optimizer=default_options.optimizer;
optim_options=rmfield(default_options,{'optimizer','blocks'});
nonlcon=default_options.nonlcon;
if isempty(nonlcon)
    nonlcon=@(x)0;
end

% get the name of the optimizer and check whether it comes from matlab, in
% which case, it is better to pass arguments as problem
if ischar(optimizer)
    optimizer_name=optimizer;
    optimizer=str2func(optimizer);
else
    optimizer_name=func2str(optimizer);
end

f0=func(x0,varargin{:});
minblk=5;
% blk_jump_prob=0.1;

npar=numel(x0);
if npar<=minblk
    blocks={1:npar};
end
nblks=numel(blocks);
if nblks==0
    if npar<=minblk
        blocks={1:npar};
        nblks=1;
    else
        nblks=ceil(npar/minblk);
        blocks=cell(nblks,1);
        randomly_rearranged=randperm(npar);
        for iblk=1:nblks
            blocks{iblk}=randomly_rearranged((iblk-1)*minblk+1:min(minblk*iblk,npar));
        end
    end
end
test=unique([blocks{:}]);
if ~isequal(test(:),(1:npar)')
    error('Not all parameters represented in the blocks')
end

Fields2Divide={'MaxIter','MaxTime','MaxFunEvals'};
for ifield=1:numel(Fields2Divide)
    v=Fields2Divide{ifield};
    optim_options.(v)=ceil(optim_options.(v)/(1+nblks));
end
% create/initialize stopping file
utils.optim.manual_stopping();
optim_options.stopping_created=true;
optim_options.funcCount=0;
optim_options.iterations=0;

blkw_obj=optim_options;
Fields2Add={'MaxIter','MaxTime','MaxFunEvals'};
for ifield=1:numel(Fields2Add)
    v=Fields2Add{ifield};
    blkw_obj.(v)=default_options.(v);
end
start_time=clock;
blkw_obj.start_time=start_time;
stopflag=utils.optim.check_convergence(blkw_obj);

matlab_style=ismember(optimizer_name,{'fmincon','fminunc','fminsearch'});
if matlab_style
    Problem=struct('objective',[],'x0',[],'Aineq',[],'bineq',[],'Aeq',[],...
        'beq',[],'lb',[],'ub',[],'nonlcon',optim_options.nonlcon,...
        'solver',optimizer_name,...
        'options',rmfield(optim_options,...
        {'nonlcon','verbose','penalty'}));
end
% add the block nonlinear constraint
optim_options.nonlcon=@block_nonlcon;
while isempty(stopflag)
    for blk=1:nblks+1
        optim_options.start_time=start_time;
        %         if rand < blk_jump_prob
        %             if rand<.5
        %                 this=random_block();
        %                 disp('optimizing a random block')
        %             else
        %                 this=1:npar;
        %                 disp('optimizing the whole vector')
        %             end
        %         else
        if blk<=nblks
            this=blocks{blk};
            disp('--------------------------------------------------------')
            disp(['optimizing block : ',int2str(blk),'/',int2str(nblks),' : ',int2str(numel(this)),' parameters'])
            disp('--------------------------------------------------------')
        else
            % now optimize the whole vector
            this=1:npar;
            disp('--------------------------------------------------------')
            disp(['optimizing block the entire vector: ',int2str(numel(this)),' parameters'])
            disp('--------------------------------------------------------')
        end
        %         end
        
        if matlab_style
            Problem.objective=@(x)block_wrapper(x,x0);
            Problem.x0=x0(this);
            Problem.lb=lb(this);
            Problem.ub=ub(this);
            Problem.options.nonlcon=@(x)block_nonlcon(x,x0);
            [x0(this),f0,exitflag,output]=optimizer(Problem);
        else
            [x0(this),f0,exitflag,output]=optimizer(@block_wrapper,x0(this),...
                lb(this),ub(this),optim_options,x0);
        end
        blkw_obj.iterations=blkw_obj.iterations+output.iterations;
        start_time=clock;
    end
    stopflag=utils.optim.check_convergence(blkw_obj);
end

    function constr=block_nonlcon(z,x)
        x(this)=z;
        constr=nonlcon(x);
    end

    function f=block_wrapper(z,x)
        x(this)=z;
        f=func(x,varargin{:});
        blkw_obj.funcCount=blkw_obj.funcCount+1;
    end

%     function blk=random_block()
%         d1=randperm(npar);
%         d2=floor(1+rand*(npar-1));
%         blk=d1(1:d2);
%     end
end

