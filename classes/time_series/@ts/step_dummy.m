function db=step_dummy(start_date,end_date,dummy_start_date)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% implements a step dummy in the time series
% inputs
%-------
% start_date: start date of the time series
% end_date: end date of the time series
% dummy_start_date: start date of the step dummy
% output
%-------
% db: a time series with zeros from the first observation to the date
% before the start of the dummy and ones from then on.

start_date=date2serial(start_date);
if ischar(dummy_start_date)
     dummy_start_date=char2serial(dummy_start_date);
elseif isnumeric(dummy_start_date)
    dummy_start_date=date2serial(dummy_start_date);
end

if numel(dummy_start_date)>1
    error('Number of dates for step dummies cannot exceed 1')
end
md=dummy_start_date;
if isempty(end_date)
    start_date=start_date:md;
else
    start_date=start_date:date2serial(end_date);
end

nobs=numel(start_date);
data=zeros(nobs,1);

loc=find(start_date==dummy_start_date);
if isempty(loc)
    error(['date ',serial2date(dummy_start_date),' could not be located in the sample'])
end
data(loc:end)=1;

db=ts(start_date,data);
end


