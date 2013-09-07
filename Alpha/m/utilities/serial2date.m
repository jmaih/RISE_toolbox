function [dat,frequency,year,period]=serial2date(s,silent)
if nargin<2
    silent=false;
end
% test=(date2serial(1990):date2serial(1995))
% serial2date(test)
% test=(date2serial('1990'):date2serial('1995'))
% serial2date(test)
% test=(date2serial('1990h1'):date2serial('1995h2'))
% serial2date(test)
% test=(date2serial('1990q1'):date2serial('1995q4'))
% serial2date(test)
% test=(date2serial('1990m1'):date2serial('1995m6'))
% serial2date(test)

[stamp,unstamp]=time_frequency_stamp();

freq=unstamp(s(1)-floor(s(1)));
year=floor(s/freq);
period=round(s-year*freq+1-stamp(freq));
period_str=num2str(period(:));
year_str=num2str(year(:));
switch freq
    case 1
        dat=cellstr(year_str);
        frequency='';
    case 2
        dat=cellstr(strcat(year_str,'H',period_str));
        frequency='H';
    case 4
        dat=cellstr(strcat(year_str,'Q',period_str));
        frequency='Q';
    case 12
        dat=cellstr(strcat(year_str,'M',period_str));
        frequency='M';
    otherwise
        if silent
            dat=[];
            frequency=[];
            return
        else
            error('wrong frequency of the time series')
        end
end
dat=cellfun(@(x)x(~isspace(x)),dat,'uniformOutput',false);