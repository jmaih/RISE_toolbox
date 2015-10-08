function y00=crs_filter_simul_init(y0,shocks,cond_shocks_id,data_y,t,p,...
    obs_id,match_first_page,nsteps)
% crs_filter_simul_initial_conditions -- initial conditions for conditional
% forecasting in the crs filters
%
% Syntax
% -------
% ::
%
%   y00=crs_filter_simul_init(y0,shocks,cond_shocks_id,data_y,t,p,...
%    obs_id,match_first_page,nsteps)
%
% Inputs
% -------
%
% - **y0** []:
%
% - **shocks** []:
%
% - **cond_shocks_id** []:
%
% - **data_y** []:
%
% - **t** []:
%
% - **p** []:
%
% - **obs_id** []:
%
% - **match_first_page** [true|{false}]:
%
% - **nsteps** []:
%
% Outputs
% --------
%
% - **y00** []:
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

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
% initial conditions for endogenous
ycond=reshape(data_y(:,t,1:nsteps),p,[]);
if match_first_page
    % all shocks can be used in the first period
    %--------------------------------------------
    shocks(:,1)=nan;
else
    % then we are not conditioning on the "true" data. In that case we
    % shift the database. We assume that we are still matching the same
    % restrictions beyond the last period
    ycond=[ycond(:,2:end),ycond(:,end)];
end
% beyond the first period, only the shocks with long reach can be used
%----------------------------------------------------------------------
shocks(cond_shocks_id,2:nsteps)=nan;

ycond=struct('data',ycond(:,:,ones(3,1)),'pos',obs_id);
econd=[shocks,zeros(exo_nbr,nsteps-1)];
econd=struct('data',econd(:,:,ones(3,1)),'pos',1:exo_nbr);
rcond=struct('data',ones(nsteps,1),'pos',nan);
y00=struct('y',y0,'ycond',ycond,'econd',econd,'rcond',rcond);
end