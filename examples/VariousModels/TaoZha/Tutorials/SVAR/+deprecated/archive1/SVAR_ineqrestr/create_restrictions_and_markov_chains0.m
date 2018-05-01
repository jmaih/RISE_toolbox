function [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains0(tpl)
% create_restrictions_and_markov_chains0 -- creates restrictions for
% the constant-parameter SVAR model
%
% ::
%
%
%   [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains0(tpl)
%
% Args:
%
%    - **tpl** [struct]: template created for SVAR objects
%
% Returns:
%    :
%
%    - **lin_restr** [cell]: two column-cell (see below). The first column
%    contains COEF objects or linear combinations of COEF objects, which are
%    themselves COEF objects.
%
%    - **nonlin_restr** [cell]: one column-cell containing inequality or
%    nonlinear restrictions
%
%    - **tpl** [struct]: non-modified template
%
% Note:
%
%    - The syntax to construct a basic COEF object is a=coef(eqtn,vname,lag)
%      - **eqtn** [integer|char]: integer or variable name
%      - **vname** [integer|char]: integer or variable name
%      - **lag** [integer]: integer or variable name
%
%    - RISE sorts the endogenous variables alphabetically and use this order
%    to tag each equation in SVAR and RFVAR models.
%
% Example:
%    coef('pi','ygap',0)
%
%    coef(2,3,0)
%
%    coef(2,'ygap',0)
%
%    coef('pi',3,0)
%
%    See also:

% syntax is coef(eqtn,vname,lag)
%-------------------------------
lin_restr={
    % first equation or "FFR" equation
    %----------------------------------
    coef('FFR','pi',1),0
    coef('FFR','pi',2),0
    coef('FFR','ygap',1),0
    coef('FFR','ygap',2),0
    coef('FFR','FFR',2),0
    % second equation or "pi" equation
    %----------------------------------
    coef('pi','FFR',0),0
    coef('pi','FFR',1),0
    coef('pi','FFR',2),0
    coef('pi','ygap',1),0
    coef('pi','ygap',2),0
    % third equation or "ygap" equation
    %-----------------------------------
    coef('ygap','FFR',1),0
    coef('ygap','FFR',2),0
    coef('ygap','pi',1),0
    coef('ygap','pi',2),0
    coef('ygap','pi',0)+coef('ygap','FFR',0),0
    };
nonlin_restr={
    % third equation or "ygap" equation
    %-----------------------------------
    'coef(ygap,FFR,0)>=0'
    };
end