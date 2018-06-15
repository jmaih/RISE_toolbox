function db=zeros(start_date,varargin)
% INTERNAL FUNCTION
%

data=zeros(varargin{:});

db=ts(start_date,data);
