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
% More About
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
% test=(date2serial('2040q1'):date2serial('2050q4'))
% serial2date(test)
% test=(date2serial('1990m1'):date2serial('1995m6'))
% serial2date(test)

% check whether it is truly serial
%-----------------------------------
if ~is_serial(s)
    
    error('input is not serial')
    
end

[stamp,unstamp]=time_frequency_stamp();

freq=unstamp(s-floor(s));

year=floor(s./freq);

period=round(s-year.*freq+1-stamp(freq));

period_str=cellstr(num2str(period(:)));

year_str=cellstr(num2str(year(:)));

dat=year_str;

frequency=cellstr(frequency2char(freq));

nonyear=freq>1;

if any(nonyear)
    
    dat(nonyear)=strcat(dat(nonyear),frequency(nonyear),period_str(nonyear));
    
end

dat=cellfun(@(x)x(~isspace(x)),dat,'uniformOutput',false);

if all(strcmp(frequency,frequency{1}))
    
    frequency=frequency{1};
    
end

dat=reshape(dat,size(s));