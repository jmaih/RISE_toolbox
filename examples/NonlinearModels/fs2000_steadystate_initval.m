function [ss,newp,retcode]=fs2000_steadystate_initval(obj,ss,pp,d,id) %#ok<INUSL>
% fs2000_steadystate_initval --  gives good start values for the steady
% state of fs2000: this is the equivalent of dynare's initvals
%
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

tmp={
    'P'			,	    2.258154387910923
    'R'			,	    1.021212121212121
    'W'			,	    4.595903784741778
    'c'			,	    0.447710752379204
    'd'			,	    0.849424911502341
    'dA'		,	    1.003004504503377
    'e'			,	    1.000000000000000
    'gp_obs'	,	    1.007971544953910
    'gy_obs'	,	    1.003004504503377
    'k'			,	    5.801216036088844
    'l'			,	    0.860424911502341
    'm'			,	    1.011000000000000
    'n'			,	    0.187215605852959
    'y'			,	    0.580765090448550
    };

retcode=0;
if nargin==1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    ss=tmp(:,1);
    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={};
else
    % no parameters to update or create in the steady state file
    %-----------------------------------------------------------
    newp=[];
    
    ys =cell2mat(tmp(:,2));
    
    % push the calculations
    %----------------------
    ss(id)=ys;
end

end
