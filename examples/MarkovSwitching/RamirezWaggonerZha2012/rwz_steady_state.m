function [y,newp,retcode]=rwz_steady_state(obj,y,p,d,id) %#ok<INUSL>
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
        y={'PAI','Y','R'};
    % flags on the calculation
    %--------------------------
    newp=struct('unique',true,'imposed',true);
else
    newp=[];
    
    PAI=1;
    Y=(p.eta-1)/p.eta;
    R=exp(p.mu_bar)/p.betta*PAI;
    
    ss=[PAI,Y,R]' ;
    
    % check the validity of the calculations
    %----------------------------------------
    if ~utils.error.valid(ss)
        retcode=1;
    else
        % push the calculations
        %----------------------
        y(id)=ss;
    end
end