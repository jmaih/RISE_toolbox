function [x0,f0,exitflag,output]=blockwise_optimization(objective,x0,lb,ub,options,varargin)

if isstruct(objective)
    x0=objective.x0;
    lb=objective.lb;
    ub=objective.ub;
    options=objective.options;
    objective=objective.objective;
end
if ischar(objective)
    objective=str2func(objective);
end

default_options=struct('MaxNodes',20,'MaxIter',1000,...
    'MaxTime',3600,'MaxFunEvals',inf,'rand_seed',100*sum(clock),...
    'penalty',1e+8,'verbose',10,'blocks',[],'optimizer',@bee_gate);
if nargin<4
    options=[];
end
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

f0=objective(x0,varargin{:});
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
manual_stopping();
optim_options.stopping_created=true;
optim_options.funcCount=0;
optim_options.iterations=0;

blkw_obj=optim_options;
Fields2Add={'MaxIter','MaxTime','MaxFunEvals'};
for ifield=1:numel(Fields2Add)
    v=Fields2Add{ifield};
    blkw_obj.(v)=default_options.(v);
end
stopflag=check_convergence(blkw_obj);
while isempty(stopflag)
    for blk=1:nblks+1
        optim_options.start_time=clock;
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
            disp(['optimizing block : ',int2str(blk),'/',int2str(nblks)])
        else
            % now optimize the whole vector
            this=1:npar;
            disp('optimizing block the entire vector')
        end
        %         end
        [x0(this),f0,exitflag,output]=optimizer(@wrapper,x0(this),lb(this),ub(this),...
            optim_options,...
            x0);
        blkw_obj.funcCount=blkw_obj.funcCount+output.funcCount;
        blkw_obj.iterations=blkw_obj.iterations+obj.iterations;
    end
    stopflag=check_convergence(blkw_obj);
end


    function f=wrapper(z,x)
        x(this)=z;
        f=objective(x,varargin{:});
    end

%     function blk=random_block()
%         d1=randperm(npar);
%         d2=floor(1+rand*(npar-1));
%         blk=d1(1:d2);
%     end
end

