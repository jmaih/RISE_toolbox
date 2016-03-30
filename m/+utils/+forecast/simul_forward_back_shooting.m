function [sims,myshocks,PAI,retcode]=simul_forward_back_shooting(T,ss,y0,state_vars_location,h,nx,Qfunc,regimes,options)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

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

% add a number of pages for the anticipation
%--------------------------------------------
myshocks=myshocks(:,:,ones(1,1+horizon));

myshocks(:,:,2:end)=0;

for istep=2:options.nsteps+1
    
    j0=istep-1;
    
    shocks=get_shocks(j0);
    
    atmp_ss=sims(:,j0)-ss{st};
    
    sims(:,istep)=ss{st}+T{st}*[atmp_ss(xlocs);0;shocks(:)];
    
    start_iter=1;
    
    if violation(sims(:,istep))
        % do a formal forward back shoot
        [tmp,shocks,start_iter]=forward_back_shoot();
        
        if retcode
            
            error(decipher(retcode))
            
        end
        
        sims(:,istep)=tmp(:,1);
        
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

    function [y1,updated_shocks,start_iter]=forward_back_shoot()
        
        is_violation=true;
        
        start=j0;
        
        shocks_=get_shocks(start,true);
        
        % these shocks are used when testing for extra violations
        these_other_shocks=zeros(size(shocks_));
            
        start_iter=0;
        
        while is_violation && start_iter<horizon+1
            
            start_iter=start_iter+1;  
            
            % set the endogenous shocks to nan over the conditioning period
            %--------------------------------------------------------------
            these_shocks=shocks_;
            
            these_shocks(is_endogenous,1:start_iter)=nan;
            
            y0_=initial_conditions_frwrd_back_shoot(these_shocks,sims(:,j0),start_iter);
            
            % forecast only for the required number of periods
            opt=utils.miscellaneous.reselect_options(options,@utils.forecast.rscond.forecast);
            opt.nsteps=start_iter;
            
            [updated_shocks,~,PAI,retcode,cfkst]=utils.forecast.rscond.forecast(model,y0_.y,...
                y0_.ycond,y0_.econd,opt,regimes(1:opt.nsteps));
            
            % remove one period of history
            y1=cfkst(:,2:start_iter+1,:);
            
            if retcode
                
                return
                
            end
            
            % The steps up until start_iter should check and so simply
            % check the next step
            
            f__=y1(:,end)-ss{st};
            
            a_expect=ss{st}+T{st}*[f__(xlocs);0;these_other_shocks(:)];
            
            is_violation=violation(a_expect);
            
        end
        
        if is_violation
            
            retcode=701;
            
        end
        
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

end