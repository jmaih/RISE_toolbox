%--- help for datetime/dateshift ---
%
% DATESHIFT Shift datetimes or generate sequences according to a calendar rule.
%    T2 = DATESHIFT(T,'start',UNIT) shifts each datetime in the array T back to
%    the beginning of the unit of time specified by UNIT. UNIT is 'year',
%    'quarter', 'month', 'week', 'day', 'hour', 'minute', or 'second'. T2 is a
%    datetime array.
% 
%    T2 = DATESHIFT(T,'end',UNIT) shifts each datetime in the array T ahead to
%    the end of the unit of time specified by UNIT. UNIT is 'year', 'quarter',
%    'month', 'week', 'day', 'hour', 'minute', or 'second'. The end of a day,
%    hour, minute, or second is also the beginning of the next one. The end of a
%    year, quarter, month, or week is the last day in that time period. T2 is a
%    datetime array.
% 
%    T2 = DATESHIFT(T,'dayofweek',DOW) returns the next occurrence of the
%    specified day of the week on or after each datetime in the array T. DOW is a
%    day of week number, or a localized day name. T2 is a datetime array.
% 
%    T2 = DATESHIFT(T,'dayofweek','weekday') returns the next occurrence of a weekday
%    on or after each datetime in the array T. DATESHIFT(T,'dayofweek','weekend')
%    returns the next occurrence of a weekend day on or after each datetime in the
%    array T.
% 
%    T2 = DATESHIFT(T,...,RULE) shifts the datetimes in the array T ahead or
%    back according to RULE. RULE is one of 'next', 'previous', or 'nearest'.
%    RULE can also be 'current' to shift to the start or end of the current
%    unit of time, or to the specified day in the current week.
% 
%    DATESHIFT treats the current day as the "next" or "previous" occurrence of
%    the specified day of week if it falls on that day of the week.
% 
%    RULE can also be an integer value or an array of integer values. For unit of
%    time, 0 corresponds to the start/end of the current unit for each datetime,
%    1 corresponds to the next, -1 to the previous, etc. For day of the week, 0
%    corresponds to the specified day in the current week for each datetime, 1
%    corresponds to the next occurrence of the specified day, -1 to the previous,
%    etc. T and RULE are the same size, or either one is a scalar.
% 
%    Examples:
% 
%       % Create an array of datetimes, and find the first and last day of the
%       % month for each one.
%          t = datetime(2013,10,30:33,10,30,0)
%          BoM = dateshift(t,'start','month')
%          EoM = dateshift(t,'end','month')
% 
%       % Create sequences of datetimes on the first and last days of each of the
%       % next six months.
%          t = datetime('now')
%          BoM = dateshift(t,'start','month',1:6)
%          EoM = dateshift(t,'end','month',1:6)
% 
%       % Create an array of datetimes, and find the first Friday on or
%       % after each one.
%          t = datetime(2013,11,20:24,10,30,0,'Format','eee, dd-MMM-yyyy HH:mm:ss')
%          fridays = dateshift(t,'dayofweek','friday')
% 
%       % Create a sequence of datetimes on the next five Fridays.
%          t = datetime('now','Format','eee, dd-MMM-yyyy HH:mm:ss')
%          fridays = dateshift(t,'dayofweek','friday',1:5)
% 
%       % Shift an array of datetimes ahead by one month, and then shift ahead to
%       % a weekday if necessary.
%          t = datetime(2013,11,20:24,10,30,0,'Format','eee, dd-MMM-yyyy HH:mm:ss')
%          monthAhead = t + calmonths(1)
%          monthAheadWeekday = dateshift(monthAhead,'dayofweek','weekday')
% 
%    See also BETWEEN, COLON.
%
%    Documentation for datetime/dateshift
%       doc datetime/dateshift
%
%    Other uses of dateshift
%
%       codistributed/dateshift    tall/dateshift
%