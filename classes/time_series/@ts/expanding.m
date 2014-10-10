function db=expanding(db,func,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


window=[];

db=ts_roll_or_expand(db,func,window,varargin{:});

end