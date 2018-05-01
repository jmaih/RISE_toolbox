function [restr,tpl]=create_restrictions_and_markov_chains0(tpl)
% create_restrictions_and_markov_chains0 -- creates linear restrictions for
% the constant-parameter SVAR model
%
% ::
%
%
%   [restr,tpl]=create_restrictions_and_markov_chains0(tpl)
%
% Args:
%
%    - **tpl** [struct]: template created for SVAR objects
%
% Returns:
%    :
%
%    - **restr** [cell]: two column-cell (see below). The first column
%    contains COEF objects or linear combinations of COEF objects, which are
%    themselves COEF objects.
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
restr={
    % first equation or "FFR" equation
    %----------------------------------
    coef(1,'pi',1),0
    coef(1,'pi',2),0
    coef(1,'ygap',1),0
    coef(1,'ygap',2),0
    coef(1,'FFR',2),0
    % second equation or "pi" equation
    %----------------------------------
    coef(2,'FFR',0),0
    coef(2,'FFR',1),0
    coef(2,'FFR',2),0
    coef(2,'ygap',1),0
    coef(2,'ygap',2),0
    % third equation or "ygap" equation
    %-----------------------------------
    coef(3,'FFR',1),0
    coef(3,'FFR',2),0
    coef(3,'pi',1),0
    coef(3,'pi',2),0
    coef(3,'pi',0)+coef(3,'FFR',0),0
    };
end