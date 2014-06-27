function [s,frequency]=date2serial(dat,silent)
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

if isnumeric(dat)
    dat=num2str(dat(:));
end
if ischar(dat)
    dat=cellstr(dat);
end
dat=upper(dat);
tmp=dat{1};
test=regexp(tmp,'\d+?(Q|M|H)\d+|\d+','tokens');
test=[test{:}];
if isempty(test)
    if any(isstrprop(tmp,'alpha'))
        if silent
            frequency=[];
            s=[];
            return
        end
        error('wrong date format')
    end
    % annual data
    freq=1;
    period=1;
    year=str2num(char(dat));
    frequency='';
else
    tmp=regexp(dat,['(?<year>\d+)',test{1},'(?<period>\d+)'],'names');
    tmp=[tmp{:}];
    year=str2num(char({tmp.year})); %#ok<*ST2NM>
    period=str2num(char({tmp.period})); %#ok<ST2NM>
    switch test{1}
        case 'H'
            freq=2;
            frequency='H';
        case 'Q'
            freq=4;
            frequency='Q';
        case 'M'
            freq=12;
            frequency='M';
    end
    if any(period<1|period>freq)||any(floor(period)-period>0)
        if silent
            frequency=[];
            s=[];
            return
        end
        error(['periods should be integers between 1 and ',sprintf('%0.0f',freq)])
    end
end
stamp=time_frequency_stamp();
dn=@(x,per,freq)x*freq+per-1+stamp(freq);

s=dn(year,period,freq);


