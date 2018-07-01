function [LogLik,Incr,retcode,Filters]=myKnownRegimesFilter(syst,y,U,z,...
    options,regs)

% myKnownRegimesFilter - filtering procedure for state-space models with
% known switching dates
%
% ::
%
%
%   [LogLik,Incr,retcode,Filters]=myKnownRegimesFilter(syst,y,U,z,options,regs)
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
%    - **regs** [vector]: history of regimes, must be of the same length as
%       **y**
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
%    - If the filter is run on a constant parameter model, then the last
%      input argument need not be specified.
%
%    - If the function is passed through a rise/dsge object, then it should
%      be called as
%      ff=filter(m,'kf_user_algo',{@myKnownRegimesFilter,regs}). That is,
%      RISE will provide all input arguments expect the last one that needs
%      to be provided by the user.
%
%    - It is the responsibility of the user to make sure that the solution
%      of the model is consistent with agents expecting no transition from
%      one regime to another. This occurs when the original model is
%      backward-looking and/or the transition matrix is diagonal.
%
%    See also: myConstantParamFilter

if ~isempty(U)
    
    error('This function still needs to be updated for non-empty U')
    
end

[T,R,Z,H,Q,sstate,init,growth]=dsge.state_space_wrapper(syst);

[~,m,h]=size(T);

n=size(y,2);

if nargin<6
    
    if h==1
        
        regs=ones(1,n);
        
    end
    
end

ybar=sstate(Z,:);

dy=ybar-sstate(Z,:);

ca=sstate;

for ireg=1:h
    
    ca(:,ireg)=(eye(m)-T(:,:,ireg))*sstate(:,ireg);
    
    ca(:,ireg)=ca(:,ireg)+growth(:,ireg);
    
end

level=options.kf_filtering_level;

nstep=options.kf_nsteps;

if numel(regs)~=n
    
    error(['the number of regimes should be equal to the number ',...
        'of observations, which is ',int2str(size(y,2))])
    
end

if any(regs<1)||any(regs>h)||~all(regs==floor(regs))
    
    error(['all regimes should be between 1 and ',int2str(h)])
    
end

init.a=init.a(:,1);

init.P=init.P(:,:,1);

f=obsw.missing_observations_kalman_filter(y,...
    T(:,:,regs),...
    R(:,:,regs),....
    Z,....
    H(:,:,regs),....
    Q(:,:,regs),init,dy,ca,level,nstep);

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
    
    myPAI=[];
    
    % assume, for filtering purposes, that the last regime will last one
    % more period
    regs_=[regs(:).',regs(end)];
    % variances (P terms) are put into cells so that format_filters
    % recognizes them quickly
    [Filters.a,Filters.P]=format_filters(n+1,permute(f.a,[1,3,2]),{f.P});
    
    Filters.PAI=myPAI;
    
    Filters.Q=repmat(eye(h),[1,1,n+1]);
    
    if level>1 % store updates
        
        [Filters.att,Filters.Ptt]=format_filters(n,...
            permute(f.att,[1,3,2]),{f.Ptt});
        
        Filters.PAItt=myPAI;
        
        if level>2 % store smoothed
            
            [Filters.atT,Filters.PtT,Filters.eta,Filters.epsilon]=format_filters(n,...
                permute(f.alpha,[1,3,2]),...
                {f.V},...
                permute(f.eta,[1,3,2]),... % smoothed shocks
                permute(f.epsilon,[1,3,2])); % smoothed measurement errors
            
            Filters.PAItT=myPAI;
            
        end
        
    end
    
end

%--------------------------
    function [varargout]=format_filters(nn,varargin)
        
        k=length(varargin);
        
        % what are the values of the unobserved when a regime is not
        % visited? we are going to assume the values are 0 instead of nan
        %------------------------------------------------------------------
        varargout=cell(1,k);
        
        is_variance=false(1,k);
        
        for iarg=1:k
            
            aa=varargin{iarg};
            
            is_variance(iarg)=iscell(aa);
            
            if is_variance(iarg)
                
                aa=aa{1};
                
                varargin{iarg}=aa;
                
            end
            
            varargout{iarg}=repmat({zeros(size(aa))},1,h);
            
        end
        
        myPAI=zeros(h,nn);
        
        for tt=1:nn
            
            myPAI(regs_(tt),tt)=1;
            
            for iarg=1:k
                
                if is_variance(iarg)
                    % should check what PP looks like when there multiple steps...
                    %--------------------------------------------------------------
                    varargout{iarg}{regs_(tt)}(:,:,tt,:)=varargin{iarg}(:,:,tt,:);
                    
                else
                    
                    varargout{iarg}{regs_(tt)}(:,:,tt)=varargin{iarg}(:,:,tt);
                    
                end
                
            end
            
        end
        
    end

end