function [y,newp,retcode]=rwz_steady_state(obj,y,p,d,id) %#ok<INUSL>
% sstate_model -- shows the way of writing a RISE steady state file
%
% ::
%
%
%   [y,newp,retcode]=sstate_model(obj,y,p,d,id)
%
% Args:
%
%    - **obj** [rise|dsge]: model object (not always needed)
%
%    - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
%    - **p** [struct]: parameter structure
%
%    - **d** [struct]: definitions
%
%    - **id** [vector]: location of the variables to calculate
%
% Returns:
%    :
%
%      CASE 1: one input argument
%
%    - **y** [cell array]: list of the variables for which the steady state
%    will be calculated within the steady state function
%
%    - **newp** [cell array]: List of the parameters calculated inside the
%    steady state function
%
%    - **retcode** []: not needed
%
%      CASE 2: More than one input argument
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
%    - If the user knows the steady state, it is always an advantage. If the
%    steady state is computed numerically, we don't know whether it is unique
%    or not. Not that it really matters but... some economists have a strong
%    aversion towards models with multiple equilibria.
%
%    - If the user does not know the solution for all the variables in the
%    steady state, it is a good idea to take a log-linear approximation for
%    the variables that potentially have a nonzero steady state. Hence the
%    user should give that information to RISE.
%
%    - One can potentially improve on the above point by explicit bounds on
%    the variables. But this is not implemented.
%
%    - An alternative that potentially avoids taking a loglinearization is to
%    to reset the values proposed by the optimizer whenever they are in a bad
%    region. It is unclear whether this always works.
%
%    - So be on the safe side, i.e. don't do like me: compute your steady
%    state analytically.
%
% Example:
%
%    See also:

retcode=0;
if nargin==1
        y={'PAI','Y','R'};
		
    newp={};
	
else
    newp=[];
    
    PAI=1;
    Y=(p.eta-1)/p.eta;
    R=exp(p.mu)/p.betta*PAI;
    
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