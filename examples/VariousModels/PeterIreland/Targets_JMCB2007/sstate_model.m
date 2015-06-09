function [y,newp,retcode]=sstate_model(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the new way of writing a RISE steady state file
%
% Syntax
% -------
% ::
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Inputs
% -------
%
% - **obj** [rise|dsge]: model object (not always needed)
%
% - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
% - **p** [struct]: parameter structure
%
% - **d** [struct]: definitions
%
% - **id** [vector]: location of the variables to calculate
%
% Outputs
% --------
%
% - **y** []: endo_nbr x 1 vector of updated steady state
%
% - **newp** [struct]: structure containing updated parameters if any
%
% - **retcode** [0|number]: return 0 if there are no problems, else return
%   any number different from 0
%
% More About
% ------------
%
% - this is new approach has three main advantages relative to the previous
%   one:
%   - The file is valid whether we have many regimes or not
%   - The user does not need to know what regime is being computed
%   - It is in sync with the steady state model
%
% Examples
% ---------
%
% See also:

retcode=0;
if nargin==1
    y={'A','THETA','Z','PAISTAR','V','PAI','Y','C','GY','GPAI','GR',...
        'LAMBDA','Q','X','R','RRPAI','E'};
    % flags on the calculation
    %--------------------------
    newp=struct('unique',true,'imposed',true,'initial_guess',false);
else
    % no parameter to update
    %-----------------------
    newp=struct();
    
    A=1;
    THETA=p.thetass;
    Z=p.zss;
    PAISTAR=1;
    V=1;
    PAI=1;
    Y=(p.thetass-1)/p.thetass*(p.zss-p.beta*p.gam)/(p.zss-p.gam);
    C=Y;
    GY=Z;
    GPAI=1;
    GR=1;
    LAMBDA=p.thetass/(p.thetass-1);
    Q=(p.zss-p.beta*p.gam)/(p.zss-p.gam);
    X=(p.thetass-1)/p.thetass;
    R=p.zss/p.beta;
    RRPAI=p.zss/p.beta;
    E=p.ess;
    
    ys=[A,THETA,Z,PAISTAR,V,PAI,Y,C,GY,GPAI,GR,LAMBDA,Q,X,R,RRPAI,E];
    ys=ys(:);
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ys)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ys;
    end
end

end