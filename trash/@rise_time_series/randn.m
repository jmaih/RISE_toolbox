function ts=randn(start_date,varargin)

data=randn(varargin{:});

ts=rise_time_series(start_date,data);

