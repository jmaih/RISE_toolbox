% DATESTR Character vector representation of date.
%   DATESTR is not recommended. Use string or char on a datetime instead.
% 
%    NOTE: While DATESTR returns text representing dates and times, it is
%    recommended that you use datetime values to represent points in time,
%    and then convert them to text using the char, string, or cellstr functions.
% 
%    S = DATESTR(V) converts one or more date vectors V to text. S is a
%    character vector or character array with M rows, where M is the number
%    of rows in V. Input V must be an M-by-6 matrix containing M full
%    (six-element) date vectors. Each element of V must be a positive
%    double-precision number.
% 
%    S = DATESTR(N) converts one or more serial date numbers N to text. S
%    has M rows, where M is the number of rows in N. Input argument N can be
%    a scalar, vector, or multidimensional array of positive
%    double-precision numbers.
% 
%    S = DATESTR(D, F) converts one or more dates into text using the format
%    F. Input argument F is a format number, character vector, or string
%    scalar that determines the format of the date character vector output.
%    D is one or more date vectors, serial date numbers, or a character
%    array or string array representing dates. Valid values for F are given
%    in Table 1, below. Input F may also contain a free-form date format
%    consisting of format tokens as shown in Table 2, below.
% 
%    Dates in a text format using 2-character years are interpreted to be
%    within the 100 years centered around the current year.
% 
%    S = DATESTR(S1, F, P) converts one or more dates into text using format
%    F. S is a character array using format F and pivot year P as the
%    starting year of the 100-year range in which a two-character year
%    resides. The default pivot year is the current year minus 50 years. F =
%    -1 uses the default format.
% 
%    S = DATESTR(...,'local') returns S in a localized format. The default
%    (which can be called with 'en_US') is US English. This argument must
%    come last in the argument sequence.
% 
%    Note:  The vectorized calling syntax can offer significant performance
%    improvement for large arrays.
% 
%    Table 1: Standard MATLAB Date format definitions
% 
%    Number           Format                   Example
%    ===========================================================================
%       0             'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17
%       1             'dd-mmm-yyyy'            01-Mar-2000
%       2             'mm/dd/yy'               03/01/00
%       3             'mmm'                    Mar
%       4             'm'                      M
%       5             'mm'                     03
%       6             'mm/dd'                  03/01
%       7             'dd'                     01
%       8             'ddd'                    Wed
%       9             'd'                      W
%      10             'yyyy'                   2000
%      11             'yy'                     00
%      12             'mmmyy'                  Mar00
%      13             'HH:MM:SS'               15:45:17
%      14             'HH:MM:SS PM'             3:45:17 PM
%      15             'HH:MM'                  15:45
%      16             'HH:MM PM'                3:45 PM
%      17             'QQ-YY'                  Q1-96
%      18             'QQ'                     Q1
%      19             'dd/mm'                  01/03
%      20             'dd/mm/yy'               01/03/00
%      21             'mmm.dd,yyyy HH:MM:SS'   Mar.01,2000 15:45:17
%      22             'mmm.dd,yyyy'            Mar.01,2000
%      23             'mm/dd/yyyy'             03/01/2000
%      24             'dd/mm/yyyy'             01/03/2000
%      25             'yy/mm/dd'               00/03/01
%      26             'yyyy/mm/dd'             2000/03/01
%      27             'QQ-YYYY'                Q1-1996
%      28             'mmmyyyy'                Mar2000
%      29 (ISO 8601)  'yyyy-mm-dd'             2000-03-01
%      30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517
%      31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17
% 
%    Table 2: Date format symbolic identifiers (Examples are in US English)
% 
%    Symbol  Interpretation of format symbol
%    ===========================================================================
%    yyyy    full year, e.g. 1990, 2000, 2002
%    yy      partial year, e.g. 90, 00, 02
%    mmmm    full name of the month, according to the calendar locale, e.g.
%            "March", "April".
%    mmm     first three letters of the month, according to the calendar
%            locale, e.g. "Mar", "Apr".
%    mm      numeric month of year, padded with leading zeros, e.g. ../03/..
%            or ../12/..
%    m       capitalized first letter of the month, according to the
%            calendar locale; for backwards compatibility.
%    dddd    full name of the weekday, according to the calendar locale, e.g.
%            "Monday", "Tuesday".
%    ddd     first three letters of the weekday, according to the calendar
%            locale, e.g. "Mon", "Tue".
%    dd      numeric day of the month, padded with leading zeros, e.g.
%            05/../.. or 20/../..
%    d       capitalized first letter of the weekday; for backwards
%            compatibility
%    QQ      quarter of year, e.g. Q2
%    HH      hour of the day, according to the time format. In case the time
%            format AM | PM is set, HH does not pad with leading zeros. In
%            case AM | PM is not set, display the hour of the day, padded
%            with leading zeros. e.g 10:20 PM, which is equivalent to 22:20;
%            9:00 AM, which is equivalent to 09:00.
%    MM      minutes of the hour, padded with leading zeros, e.g. 10:15,
%            10:05, 10:05 AM.
%    SS      second of the minute, padded with leading zeros, e.g. 10:15:30,
%            10:05:30, 10:05:30 AM.
%    FFF     milliseconds field, padded with leading zeros, e.g.
%            10:15:30.015.
%    PM      AM or PM
% 
%    Examples:
% 	DATESTR(now) returns '24-Jan-2003 11:58:15' for that particular date,
% 	on an US English locale DATESTR(now,2) returns 01/24/03, the same as
% 	for DATESTR(now,'mm/dd/yy') DATESTR(now,'dd.mm.yyyy') returns
% 	24.01.2003 To convert a non-standard date form into a standard MATLAB
% 	dateform, first convert the non-standard date form to a date number,
% 	using DATENUM, for example,
% 	DATESTR(DATENUM('24.01.2003','dd.mm.yyyy'),2) returns 01/24/03.
% 
% 	See also DATETIME, DATETIME/STRING, DATETICK.
%
%    Documentation for datestr
%       doc datestr
%
%    Other uses of datestr
%
%       calendarDuration/datestr    rise_dates.dates/datestr
%       datetime/datestr            tall/datestr
%