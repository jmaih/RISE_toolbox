function Results=simple_rwmh_sampler(TargetFun,x0,lb,ub,Sigma_lower_chol,...
    nblks,nsimPerBlck,varargin)

Results=struct();

f0=TargetFun(x0,varargin{:});

npar=numel(x0);

block_iter=0;

last_draw=x0(:,ones(1,nblks));

last_posterior=f0*ones(1,nblks);

nruns=nsimPerBlck*ones(1,nblks);

Results.pop=struct('x',{},'f',{});

Results.AcceptanceRate=zeros(1,nblks);

for blk = 1:nblks
    
    block_iter=block_iter+1;
    
    accepted_draws_this_chain = 0;
    
    feval_this_chain = 0;
    
    draw_iter = 0;
    
    msg=sprintf('RWMH block %0.0f/%0.0f',blk,nblks);
    
    h = waitbar(draw_iter,msg);
    
    while draw_iter < nruns(blk)
        
        draw_iter=draw_iter+1;
        
        [par, logpost, accepted] = draw_one(TargetFun, last_draw(:,blk), ...
            last_posterior(blk), Sigma_lower_chol,lb,ub,npar,...
            varargin{:});
        
        last_draw(:,blk) = par;
        
        last_posterior(blk) = logpost;
        
        % store results
        %--------------
        store_results()
        
        feval_this_chain = feval_this_chain + 1;
        
        accepted_draws_this_chain = accepted_draws_this_chain + accepted;
        
        accept_rate=100*accepted_draws_this_chain/draw_iter;
        
        msg=sprintf('RWMH block %0.0f/%0.0f, accept rate %0.2f',blk,nblks,...
            accept_rate);
        
        waitbar(draw_iter/nruns(blk),h,msg)
        
    end
    
    close(h)
    
    Results.AcceptanceRate(blk) = accepted_draws_this_chain/draw_iter;
    
end

    function store_results()
        
        Results.pop(blk,draw_iter).x=par;
        
        Results.pop(blk,draw_iter).f=logpost;
        
    end

end

function  [par, logpost, accepted] = draw_one(...
    TargetFun,last_draw,last_posterior,Sigma_lower_chol,lb,ub,n,varargin)

%  case 'random_walk_metropolis_hastings'

ProposalFun = @rand_multivariate_normal;

par = ProposalFun(last_draw);

if all( par(:) > lb ) && all( par(:) < ub )
    
    try
        
        logpost = TargetFun(par(:),varargin{:});
        
    catch
        
        logpost = inf;
        
    end
    
else
    
    logpost = inf;
    
end

accepted = false;

if logpost<inf
    
    r = -(logpost-last_posterior);
    
    log_rand=log(rand);
    
    accepted =  log_rand < r;
    
end

if ~accepted
    
    par = last_draw;
    
    logpost = last_posterior;
    
end

    function draw = rand_multivariate_normal(Mean)
        
        draw = Mean + Sigma_lower_chol*randn(n,1);
        
    end

end
