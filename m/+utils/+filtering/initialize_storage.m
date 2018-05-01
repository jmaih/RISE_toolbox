function [Filters,K_store,iF_store,v_store]=initialize_storage(a,P,PAI,Q,p0,exo_nbr,horizon,nsteps,smpl,store_filters)
% initialize_storage - initialize the storage of filtering
%
% ::
%
%
%   [Filters,K_store,iF_store,v_store]=initialize_storage(a,P,PAI,Q,p0,...
%   exo_nbr,horizon,m,nsteps,smpl,h,store_filters)
%
% Args:
%
%    - **a** [cell]: initial conditions of the filter for each regime
%
%    - **P** [cell]: initial covariance of the filter for each regime
%
%    - **PAI** [vector]: initial probability distribution of regimes
%
%    - **p0** [scalar]: number of observables
%
%    - **exo_nbr** [scalar]: number of exogenous
%
%    - **horizon** [{1}|scalar]: number of anticipated steps + 1
%
%    - **nsteps** [scalar]: number of forecast steps
%
%    - **smpl** [scalar]: number of observations
%
%    - **store_filters** [0|1|2|3]: 0 (no storage), 1(predicted only),
%     2(predicted and updated), 3(predicted, updated and smoothed)
%
% Returns:
%    :
%
%    - **Filters** [struct]: structure with different fields
%
%    - **K_store** [cell]: Place holder for Kalman gains
%
%    - **iF_store** [cell]: Place holder for inverses of covariance matrices
%    of forecast errors
%
%    - **v_store** [cell]: Place holder for forecast errors
%
% Note:
%
% Example:
%
%    See also:

Filters=struct();
K_store=[];
iF_store=[];
v_store=[];
if store_filters>0 % store filters
    h=numel(a);
    m=size(a{1},1);
    Filters.a=repmat({zeros(m,nsteps,smpl+1)},1,h);
    Filters.P=repmat({zeros(m,m,smpl+1)},1,h);
    for state=1:h
        Filters.a{state}(:,1,1)=a{state};
        Filters.P{state}(:,:,1)=P{state};
    end
    Filters.PAI=zeros(h,smpl+1);
    Filters.PAI(:,1)=PAI;
    for istep=2:nsteps
        % in steady state, we remain at the steady state
        %------------------------------------------------
        for state=1:h
            Filters.a{state}(:,istep,1)=Filters.a{state}(:,istep-1,1);
        end
    end
    Filters.Q=zeros(h,h,smpl+1);
    Filters.Q(:,:,1)=Q;
    if store_filters>1 % store updates
        Filters.att=repmat({zeros(m,1,smpl)},1,h);
        Filters.Ptt=repmat({zeros(m,m,smpl)},1,h);
        Filters.PAItt=zeros(h,smpl);
        if store_filters>2 % store smoothed
            K_store=repmat({zeros(m,p0,smpl)},1,h);
            iF_store=repmat({zeros(p0,p0,smpl)},1,h);
            v_store=repmat({zeros(p0,smpl)},1,h);
            Filters.atT=repmat({zeros(m,1,smpl)},1,h);
            Filters.PtT=repmat({zeros(m,m,smpl)},1,h);
            Filters.eta=repmat({zeros(exo_nbr*horizon,1,smpl)},1,h); % smoothed shocks
            Filters.epsilon=repmat({zeros(p0,1,smpl)},1,h); % smoothed measurement errors
            Filters.PAItT=zeros(h,smpl);
        end
    end
end
end
