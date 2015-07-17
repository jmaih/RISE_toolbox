function db=zeros(start_date,varargin)
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


data=zeros(varargin{:});

db=ts(start_date,data);

