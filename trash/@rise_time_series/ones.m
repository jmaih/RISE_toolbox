function ts=ones(start_date,varargin)

data=ones(varargin{:});

ts=rise_time_series(start_date,data);

