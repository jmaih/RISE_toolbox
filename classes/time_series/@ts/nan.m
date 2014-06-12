function db=nan(start_date,varargin)

data=nan(varargin{:});

db=ts(start_date,data,[],[],true);

