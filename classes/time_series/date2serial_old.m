function [s,frequency]=date2serial(dat,silent)
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
if isempty(dat)
    s=[];
    frequency=[];
    return
elseif is_serial(dat)
    dat=serial2date(dat);
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

% check whether it is already serial
%-----------------------------------
[freq,frequency]=serial2frequency(dat);
if ~isempty(freq)
    s=dat;
else
    if isnumeric(dat)
        dat=num2str(dat(:));
    end
    if ischar(dat)
        dat=cellstr(dat);
    end
    dat=upper(dat);
    tmp=dat{1};
    fmap=frequency_map();
    WMQH=strrep(cell2mat(strcat(fmap.strings','|')),'||','');
    test=regexp(tmp,['\d+?(',WMQH,')\d+|\d+'],'tokens');
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
        period=str2num(char({tmp.period}));
        frequency=test{1};
        freq=frequency2num(frequency);
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
end

end

