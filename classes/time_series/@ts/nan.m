function db=nan(start_date,varargin)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


data=nan(varargin{:});

db=ts(start_date,data,[],[],true);

