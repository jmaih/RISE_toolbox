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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


data=zeros(varargin{:});

db=ts(start_date,data);

