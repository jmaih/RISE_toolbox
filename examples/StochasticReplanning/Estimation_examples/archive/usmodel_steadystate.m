function [y,newp,retcode]=usmodel_steadystate(obj,y,p,d,id) %#ok<INUSL>
% steady_state_file -- shows the new way of writing a RISE steady state
% file
%
% Syntax
% -------
% ::
%
%   [y,newp,retcode]=steady_state_file(obj,y,p,d,id)
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
    y={'dy','dc','dinve','dw','pinfobs','robs','labobs'};
    % flags on the calculation
    %--------------------------
    newp=struct('unique',false,'imposed',false,'initial_guess',true);
else
    % if some parameters are computed in the steady state, they have to be
    % returned in a structure or in a cell with two columns
    %----------------------------------------------------------------------
    newp=[];
    
    % In the SW model, one of the steady state is endogenous...
    cpie=1+p.constepinf/100;
    cgamma=1+p.ctrend/100 ;
    cbeta=1/(1+p.constebeta/100);
    cr=cpie/(cbeta*cgamma^(-p.csigma));
    conster=(cr-1)*100;
    
    
    dy=ctrend;
    dc=ctrend;
    dinve=ctrend;
    dw=ctrend;
    pinfobs = constepinf;
    robs =conster;
    labobs =constelab;
    
    y =[dy,dc,dinve,dw,pinfobs,robs,labobs]';
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(y)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ys;
    end
end

