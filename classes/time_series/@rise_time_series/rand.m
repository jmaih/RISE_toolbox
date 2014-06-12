function ts=rand(start_date,varargin)

data=rand(varargin{:});

ts=rise_time_series(start_date,data);

