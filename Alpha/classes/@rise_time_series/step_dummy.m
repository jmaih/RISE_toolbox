function ts=step_dummy(start_date,end_date,dummy_start_date)
% implements a step dummy in the time series
% inputs
%-------
% start_date: start date of the time series
% end_date: end date of the time series
% dummy_start_date: start date of the step dummy
% output
%-------
% ts: a time series with zeros from the first observation to the date
% before the start of the dummy and ones from then on.

if ischar(start_date)
    start_date={start_date};
end
if ischar(dummy_start_date)
    dummy_start_date={dummy_start_date};
end
dummy_start_date=rise_date(dummy_start_date);
if numel(dummy_start_date)>1
    error('multiple start dates for step dummies not allowed')
end

ts=step_dummy(rise_date(start_date),rise_date(end_date),rise_date(dummy_start_date):rise_date(end_date));

