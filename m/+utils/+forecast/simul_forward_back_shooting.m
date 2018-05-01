function [sims,myshocks,PAI,retcode]=simul_forward_back_shooting(T,ss,y0,state_vars_location,h,nx,Qfunc,regimes,options)
% H1 line
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
%      At present the restrictions are enforced only in the violating periods,
%      which implies that the contemporaneous shocks are altered. But one can
%      think of a scheme whereby only future shocks are changed, keeping the
%      original/actual shocks unchanged.
%
% Example:
%
%    See also:

cutoff=-1e-10;

PAI=1;

if h>1
    
    error('Forward-back shooting only available for constant-parameter models')
    
end

st=1;

xlocs=state_vars_location;

if options.simul_anticipate_zero
    
    error('simul_anticipate_zero inconsistent with forward-back shooting')
    
end

retcode=0;

model=struct('T',{T},'sstate',{ss},'state_cols',state_vars_location,...
    'Qfunc',Qfunc,'k',options.k_future,'nshocks',nx);

horizon=options.k_future;

if horizon==0
    
    error('Must anticipate the future for forward-back shooting')
    
end

sims=y0.y(:,ones(options.nsteps+1,1));

violation=@(x)any(options.sep_compl(x)<cutoff);

threshold_lookup();

% straigthen the shocks (they may be in a wrong order)...
%--------------------------------------------------------
myshocks=nan(nx,size(y0.econd.data,2));

myshocks(y0.econd.pos,:)=y0.econd.data(:,:,1);

is_endogenous=options.is_endogenous_exo_vars;

if ~any(is_endogenous)
    
    is_endogenous=all(isnan(myshocks),2);
    
    if ~any(is_endogenous)
        
        error(['Forward-back shooting requires the presence of ',...
            'exogenous variables to endogenize'])
        
    end

end

simul_fbs_horizon=options.simul_fbs_horizon;

if ~ismember(simul_fbs_horizon,(0:horizon))
    
    error(['simul_fbs_horizon must be an integer in [0,',int2str(horizon),']'])
    
end

howmany=simul_fbs_horizon+1;

% add a number of pages for the anticipation
%--------------------------------------------
myshocks=myshocks(:,:,ones(1,1+horizon));

myshocks(:,:,2:end)=0;

for istep=2:options.nsteps+1
    
    j0=istep-1;
    
    shocks=get_shocks(j0);
    
    [y1,isviol,icount]=do_multiple_steps(sims(:,j0),shocks,howmany);
    
    if isviol
        % do a formal forward back shoot
        [tmp,shocks,start_iter]=forward_back_shoot(icount,j0);
        
        if retcode
            
            error(decipher(retcode))
            
        end
        
        sims(:,istep)=tmp(:,1);
        
    else
        
        sims(:,istep)=y1;
        
        start_iter=1;
        
    end
    
    % replace the updated shocks
    myshocks(is_endogenous,j0,1:start_iter)=shocks(is_endogenous,1:start_iter);
    
end

myshocks=myshocks(:,1:options.nsteps,:);

% Keep history in, it will be removed by the caller
%--------------------------------------------------
sims=sims(:,1:end);

    function shocks=get_shocks(j0,keep_nan)
        
        if nargin<2
            
            keep_nan=false;
            
        end
        
        limit=min(j0+options.k_future,size(myshocks,2));
        
        shocks=myshocks(:,j0:limit);
        
        if ~keep_nan
            
            shocks(isnan(shocks))=0;
            
        end
        
        missing=options.k_future+1-size(shocks,2);
        
        if missing>0
            
            shocks=[shocks,zeros(nx,missing)]; 
            
        elseif missing<0
            
            error('something quite not right')
            
        end
        
        % set all the shocks beyond the first period to 0
        %-------------------------------------------------
        shocks(:,2:end)=0;
        
    end

    function [y1,updated_shocks,start_iter]=forward_back_shoot(icount,j0)
        
        % icount: first period of violation
        % j0: date of first shock
        
        is_violation=true;
        
        shocks_=get_shocks(j0,true);
        
        % these shocks are used when testing for extra violations
        these_other_shocks=zeros(size(shocks_));
            
        start_iter=icount-1;
        
        first_is_viol=icount==1;
        
        first_was_viol=first_is_viol;
        
        while is_violation && start_iter<horizon+1
            
            start_iter=start_iter+1;  
            
            % set the endogenous shocks to nan over the conditioning period
            %--------------------------------------------------------------
            these_shocks=shocks_;
            
            if first_is_viol
                
                istart=1;
                
            else
                
                istart=max(1,icount);
                
            end
            
            these_shocks(is_endogenous,istart:start_iter)=nan;
            
            y0_=initial_conditions_frwrd_back_shoot(these_shocks,sims(:,j0),start_iter);
            
            % update the initial conditions
            %------------------------------
            if istart>1
                % do not constrain the first period since the shock in that
                % period is given
                %-------------------------------------------------------
                y0_.ycond.data(:,1:icount-1,:)=nan;
                
            end
            
            % forecast only for the required number of periods
            opt=utils.miscellaneous.reselect_options(options,@utils.forecast.rscond.forecast);
            opt.nsteps=start_iter;
            
            [updated_shocks,~,PAI,retcode,cfkst]=utils.forecast.rscond.forecast(model,y0_.y,...
                y0_.ycond,y0_.econd,opt,regimes(1:opt.nsteps));
            
            % remove one period of history
            %------------------------------
            y1=cfkst(:,2:start_iter+1,:);
            
            if retcode
                
                return
                
            end
            
            % recheck all the forecast starting from the first since it is
            % affected by future shocks
            %--------------------------------------------------------------
            
            for icol=1:size(y1,2)
                
                is_violation=violation(y1(:,icol));
                
                if is_violation
                    
                    if icol==1 && ~first_was_viol
                        
                        first_is_viol=true;
                        
                    end
                    
                    icount=icol;
                    
                    start_iter=start_iter-1;
                    
                    break
                    
                end
                
            end
            
            % Then check an additional step
            %--------------------------------
            if ~is_violation
                
                [~,is_violation]=one_step(y1(:,end),these_other_shocks);
                
            end

        end
        
        if is_violation
            
            retcode=701;
            
        end
        
    end

    function [y1,isviol,icount]=do_multiple_steps(y0,shocks,howmany)
        % Do "howmany" steps while there is no violation
        %-----------------------------------------------
        
        if nargin<3
            
            howmany=1;
            
        end
        
        % in the first step, use the original shocks
        %-------------------------------------------
        [y1,isviol]=one_step(y0,shocks);
        
        % use zero shocks in the subsequent ones
        %----------------------------------------
        shocks=0*shocks;
        
        y1=y1(:,ones(howmany,1));
        
        icount=1;
        
        while ~isviol && icount<howmany
            
            icount=icount+1;
            
            [y1(:,icount),isviol]=one_step(y1(:,icount-1),shocks);
            
        end
        
        y1=y1(:,1);
        
    end

    function [y1,isviol]=one_step(y0,shocks)
        
        y0=y0-ss{st};
        
        y1=ss{st}+T{st}*[y0(xlocs);0;shocks(:)];
        
        isviol=violation(y1);
        
    end

    function y0_=initial_conditions_frwrd_back_shoot(the_shocks,a_start,start_iter)
        % compute a conditional forecast
        
        ycond=y0.ycond.data(:,1:start_iter,1);
        
        ycond=struct('data',ycond(:,:,ones(3,1)),'pos',y0.ycond.pos);
        
        econd=[the_shocks,zeros(nx,start_iter-1)];
        
        econd=struct('data',econd(:,:,ones(3,1)),'pos',1:nx);
        
        rcond=struct('data',ones(start_iter,1),'pos',nan);
        
        y0_=struct('y',a_start,'ycond',ycond,'econd',econd,'rcond',rcond);
        
    end

    function threshold_lookup()
        
        if ~isempty(y0.ycond.data)
            
            return
            
        end
                
        n=numel(y0.y);
        
        z0=zeros(n,1);
        
        pos=y0.ycond.pos;
        
        b0=zeros(numel(pos),1);
        
        fmsOptions=struct('Display','none');
        
        b=fminsearch(@lookup,b0,fmsOptions);
        
        y0.ycond.data=b(:,ones(1,horizon),ones(1,3));
        
        function r=lookup(b)
            
            z=z0;
            
            z(pos)=b;
            
            r=options.sep_compl(z);
            
            r=abs(r)+violation(z);
            
        end
        
    end

end