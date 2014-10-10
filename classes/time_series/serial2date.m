function [dat,frequency,year,period]=serial2date(s,silent)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

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

% check whether it is truly serial
%-----------------------------------
if ~is_serial(s)
    error('input is not serial')
end

[stamp,unstamp]=time_frequency_stamp();

freq=unstamp(s(1)-floor(s(1)));
year=floor(s/freq);
period=round(s-year*freq+1-stamp(freq));
period_str=num2str(period(:));
year_str=num2str(year(:));
dat=year_str;
fmap=frequency_map();
if ~ismember(freq,fmap.code)
    if silent
        dat=[];
        frequency=[];
        return
    else
        error('wrong frequency of the time series')
    end
end
frequency=frequency2char(freq);
if freq>1
    dat=strcat(dat,frequency,period_str);
end
dat=cellstr(dat);
dat=cellfun(@(x)x(~isspace(x)),dat,'uniformOutput',false);