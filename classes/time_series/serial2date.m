function [dat,frequency,year,period]=serial2date(s,silent)
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

dec=serial2dec(s,true);

year=[dec.year];

period=[dec.period];

period_str=cellstr(int2str(period(:)));

year_str=cellstr(int2str(year(:)));

dat=year_str;

frequency=frequency2char(dec(1).freq);

nonyear=dec(1).freq>1;

if nonyear
    
    dat=strcat(dat,frequency,period_str);
    
end

dat=cellfun(@(x)x(~isspace(x)),dat,'uniformOutput',false);

dat=reshape(dat,size(s));