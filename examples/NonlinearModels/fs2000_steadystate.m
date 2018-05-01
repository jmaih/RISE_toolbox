function [ss,newp,retcode]=fs2000_steadystate(obj,ss,pp,d,id) %#ok<INUSL>
% fs2000_steadystate --  computes the steady state of fs2000 analytically
%
% ::
%
%
%   [y,newp,retcode]=fs2000_steadystate(obj,y,p,d,id)
%
% Args:
%
%    - **obj** [rise|dsge]: model object (not always needed)
%
%    - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
%    - **pp** [struct]: parameter structure
%
%    - **d** [struct]: definitions
%
%    - **id** [vector]: location of the variables to calculate
%
% Returns:
%    :
%
%    - **y** []: endo_nbr x 1 vector of updated steady state
%
%    - **newp** [struct]: structure containing updated parameters if any
%
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0
%
% Note:
%
%    - this is new approach has three main advantages relative to the previous
%      one:
%      - The file is valid whether we have many regimes or not
%      - The user does not need to know what regime is being computed
%      - It is in sync with the steady state model
%
% Example:
%
%    See also:

retcode=0;
if nargin==1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    ss={'m','P','c','e','W','R','k','d','n','l','gy_obs','gp_obs','y','dA'};
    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={};
else
    % no parameters to update or create in the steady state file
    %-----------------------------------------------------------
    newp=[];
    
    % computation of the steady state
    %--------------------------------
    dA = exp(pp.gam);
    gst = 1/dA;
    m = pp.mst;
    
    khst = ( (1-gst*pp.bet*(1-pp.del)) / (pp.alp*gst^pp.alp*pp.bet) )^(1/(pp.alp-1));
    xist = ( ((khst*gst)^pp.alp - (1-gst*(1-pp.del))*khst)/pp.mst )^(-1);
    nust = pp.psi*pp.mst^2/( (1-pp.alp)*(1-pp.psi)*pp.bet*gst^pp.alp*khst^pp.alp );
    n  = xist/(nust+xist);
    P  = xist + nust;
    k  = khst*n;
    
    l  = pp.psi*pp.mst*n/( (1-pp.psi)*(1-n) );
    c  = pp.mst/P;
    d  = l - pp.mst + 1;
    y  = k^pp.alp*n^(1-pp.alp)*gst^pp.alp;
    R  = pp.mst/pp.bet;
    W  = l/n;
    %   ist  = y-c;
    %   q  = 1 - d;
    
    e = 1;
    
    gp_obs = log(m/dA);
    gy_obs = log(dA);
    
    ys =[m P c e W R k d n l gy_obs gp_obs y dA]';
    
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ys)
        retcode=1;
    else
        % push the calculations
        %----------------------
        ss(id)=ys;
    end
end

end
