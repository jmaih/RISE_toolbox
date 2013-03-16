function ts=zeros(start_date,varargin)

data=zeros(varargin{:});

ts=rise_time_series(start_date,data);

