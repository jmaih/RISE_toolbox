function [LogLik,Incr,retcode,Filters]=myConstantParamFilter(syst,y,U,z,options)
% myConstantParamFilter - filtering procedure for state-space models with
% constant parameters
%
% ::
%
%
%   [LogLik,Incr,retcode,Filters]=myConstantParamFilter(syst,y,U,z,options)
%
% Args:
%
%    - **syst** [struct]: structure provided by dsge.filter
%
%    - **y** [matrix]: matrix of data provided by dsge.filter
%
%    - **U** [matrix]: matrix of trends provided by dsge.filter
%
%    - **z** [matrix]: matrix of deterministic terms provided by dsge.filter
%
%    - **options** [struct]: options provided by dsge.filter
%
% Returns:
%    :
%
%    - **LogLik** [numeric]: value of the log likelihood
%
%    - **Incr** [vector]: contributions to the likelihood in each period
%
%    - **retcode** [numeric]: flag equal to 0 if there is no problem  
%
%    - **Filters** [struct]: structure containing all the filtering
%      information
%
% Note:
%
%    - If the function is passed through a rise/dsge object, then it should
%      be called as ff=filter(m,'kf_user_algo',@myConstantParamFilter).
%
%    See also: myKnownRegimesFilter

if ~isempty(U)
    
    error('This function still needs to be updated for non-empty U')
    
end

[T,R,Z,H,Q,sstate,init,growth]=dsge.state_space_wrapper(syst);

ybar=sstate(Z);

dy=ybar-sstate(Z);

m=size(T,1);

ca=(eye(m)-T)*sstate;

if any(growth~=0)
    
    ca = ca+growth;
    
end

level=options.kf_filtering_level;

nstep=options.kf_nsteps;

f=obsw.missing_observations_kalman_filter(y,T,R,Z,H,Q,init,dy,ca,level,nstep);

retcode=f.retcode;

Incr=f.incr;

% override the log lik if there is a presample...

if options.kf_presample
    
    LogLik=sum(Incr(options.kf_presample+1:end));
    
else
    
    LogLik=f.log_lik;
    
end

Filters=struct();

%--------------------------
if level>0 % store filters
    
    n=size(y,2);
    
    Filters.a={permute(f.a,[1,3,2])};
    
    Filters.P={f.P};
    
    Filters.PAI=ones(1,n+1);
    
    Filters.Q=ones(1,1,n+1);
    
    if level>1 % store updates
        
        Filters.att={permute(f.att,[1,3,2])};
        
        Filters.Ptt={f.Ptt};
        
        Filters.PAItt=ones(1,n);
        
        if level>2 % store smoothed
            
            Filters.atT={permute(f.alpha,[1,3,2])};
            
            Filters.PtT={f.V};
            
            Filters.eta={permute(f.eta,[1,3,2])}; % smoothed shocks
            
            Filters.epsilon={permute(f.epsilon,[1,3,2])}; % smoothed measurement errors
            
            Filters.PAItT=ones(1,n);
            
        end
        
    end
    
end

%--------------------------

end