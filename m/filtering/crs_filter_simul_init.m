function y00=crs_filter_simul_init(y0,shocks,cond_shocks_id,data_y,t,p,...
obs_id,match_first_page,nsteps)
% crs_filter_simul_initial_conditions -- initial conditions for conditional
% forecasting in the crs filters
%
% ::
%
%
%   y00=crs_filter_simul_init(y0,shocks,cond_shocks_id,data_y,t,p,...
%    obs_id,match_first_page,nsteps)
%
% Args:
%
%    - **y0** []:
%
%    - **shocks** []:
%
%    - **cond_shocks_id** []:
%
%    - **data_y** []:
%
%    - **t** []:
%
%    - **p** []:
%
%    - **obs_id** []:
%
%    - **match_first_page** [true|{false}]:
%
%    - **nsteps** []:
%
% Returns:
%    :
%
%    - **y00** []:
%
% Note:
%
% Example:
%
%    See also:

if nargin<9
    
    nsteps=[];
    
    if nargin<8
        
        match_first_page=false;
        
    end
    
end 

[exo_nbr,horizon]=size(shocks);

if isempty(nsteps)
    
    nsteps=horizon;
    
end

if match_first_page
    % initial conditions for endogenous. Skip first page if not matched.
    ycond=reshape(data_y(:,t,1:nsteps),p,[]);
    % all shocks can be used in the first period
    %--------------------------------------------
    shocks(:,1)=nan;
    
else
    % then we are not conditioning on the "true" data. We assume that we
    % are still matching the same restrictions throughout
    ycond=reshape(data_y(:,t,2),p,[]);
    
    ycond=ycond(:,ones(1,nsteps));
    
end
% beyond the first period, only the shocks with long reach can be used
%----------------------------------------------------------------------
shocks(cond_shocks_id,1:nsteps)=nan;

ycond=struct('data',ycond(:,:,ones(3,1)),'pos',obs_id);

econd=[shocks,zeros(exo_nbr,nsteps-1)];

econd=struct('data',econd(:,:,ones(3,1)),'pos',1:exo_nbr);

rcond=struct('data',ones(nsteps,1),'pos',nan);

y00=struct('y',y0,'ycond',ycond,'econd',econd,'rcond',rcond);

end