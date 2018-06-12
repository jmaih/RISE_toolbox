function db=zeros(start_date,varargin)

data=zeros(varargin{:});

db=ts(start_date,data);
