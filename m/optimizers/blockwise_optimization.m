function [x0,f0]=blockwise_optimization(objective,x0,lb,ub,blocks,optimfunc,varargin)

if nargin<6
    optimfunc=[];
    if nargin<5
        blocks=[];
    end
end


f0=objective(x0,varargin{:});
minblk=15;
blk_jump_prob=0.1;

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
if isempty(optimfunc)
    optimfunc=@bee;
end

% create/initialize stopping file
manual_stopping();
% obj=bee(Objective,x0,f0,lb,ub,options,varargin)
optim_options=struct('start_time',clock,'colony_size',20,'max_iter',100,...
    'max_time',3600,'max_fcount',inf,'rand_seed',100*sum(clock),...
    'penalty',1e+8,'verbose',10,'stopping_created',true,'fcount',0,'iter',0);
blkw_obj=optim_options;
blkw_obj.max_fcount=10000;
stopflag=check_convergence(blkw_obj);
while isempty(stopflag)
    for blk=1:nblks
        if rand < blk_jump_prob
            if rand<.5
                this=random_block();
                disp('optimizing a random block')
            else
                this=1:npar;
                disp('optimizing the whole vector')
            end
        else
            this=blocks{blk};
        end
        disp(['optimizing block : ',int2str(blk),'/',int2str(nblks)])
        obj=optimfunc(@wrapper,x0(this),f0,lb(this),ub(this),...
            optim_options,...
            x0);
        x0(this)=obj.best;
        f0=obj.best_fval;
        optim_options.fcount=obj.fcount;
        blkw_obj.fcount=optim_options.fcount;
        blkw_obj.iter=blkw_obj.iter+obj.iter;
    end
    stopflag=check_convergence(blkw_obj);
end


    function f=wrapper(z,x)
        x(this)=z;
        f=objective(x,varargin{:});
    end

    function blk=random_block()
        d1=randperm(npar);
        d2=floor(1+rand*(npar-1));
        blk=d1(1:d2);
    end
end

