function db=ones(start_date,varargin)

data=ones(varargin{:});

db=ts(start_date,data);

