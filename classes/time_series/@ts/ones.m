function db=ones(start_date,varargin)
% ones overloads ones for ts objects
%
% Syntax
% -------
% ::
%   db=ts.ones(start_date,varargin)
%
% Inputs
% -------
%
% - **start_date** : [numeric|char]: a valid time series (ts) date
%
% - varargin : [numeric]: arguments to matlab's **ones** function.
%
% Outputs
% --------
%
% - **db** : [ts]: a time series
%
% Description
% ------------
%
% - this is a static method and so it has to be called with the **ts.**
%   prefix
%
% - ts.ones does not allow more than 3 dimensions
%
% Examples
% ---------
%
%   db=ts.ones(1990,10,1)
%   db=ts.ones('1990',10,3)
%   db=ts.ones('1990Q3',10,5,100)
%
% See also: 


data=ones(varargin{:});

db=ts(start_date,data);

