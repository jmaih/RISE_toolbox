%  SERIAL2DATE Convert serial dates to date representations.
% 
%  [dat, frequency] = serial2date(s) converts the input serial dates s to
%  date representations. The output 'dat' is an array or cell array of
%  date strings corresponding to the input serial dates. The 'frequency'
%  output indicates the frequency of the dates.
% 
%  Inputs:
%  - s: Serial dates. It can be either an array of doubles representing
%  undated or non-dated values, or an array of objects of the
%  'rise_dates.dates' class representing dated values.
% 
%  Outputs:
%  - dat: Array or cell array of date strings corresponding to the input
%  serial dates. The size of 'dat' is the same as the size of 's'.
%  - frequency: Frequency of the dates. It can take the following values:
%  'U' (undated/non-dated), 'D' (daily), 'W' (weekly), 'M' (monthly),
%  'Q' (quarterly), 'H' (half-yearly), 'Y' (yearly), or 'A' (annual).
% 
%  Examples:
%  - Undated/Non-dated:
%  test = (1990:1995);
%  serial2date(test)
%  - Daily dates:
%  test = (date2serial('1/22/1990'):date2serial('03/25/1995'));
%  serial2date(test)
%  - Weekly dates:
%  test = (date2serial('1990w1'):date2serial('1995w16'));
%  serial2date(test)
%  - Monthly dates:
%  test = (date2serial('1990m1'):date2serial('1995m6'));
%  serial2date(test)
%  - Quarterly dates:
%  test = (date2serial('2040q1'):date2serial('2050q4'));
%  serial2date(test)
%  - Half-yearly dates:
%  test = (date2serial('1990h1'):date2serial('1995h2'));
%  serial2date(test)
%  - Annual or yearly dates:
%  test = (date2serial('1990'):date2serial('1995'));
%  serial2date(test)
%