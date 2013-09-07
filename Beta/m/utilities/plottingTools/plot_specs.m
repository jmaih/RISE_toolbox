function pp=plot_specs(TimeInfo,nticks)

% see also plot_real_time

if nargin<2
    nticks=5;
end
pp=struct();
pp.xdatenums=[TimeInfo.date_number];
pp.nticks=min(nticks,numel(pp.xdatenums));
pp.tickLocs=linspace(pp.xdatenums(1),pp.xdatenums(end),pp.nticks);
for ii=1:pp.nticks
    pp.tickLocs(ii)=find_nearest(pp.xdatenums,pp.tickLocs(ii));
end
switch TimeInfo(1).freq
    case ''
        date_format='yy';
    case {'H','Q'}
        date_format='QQ-YY';
    case {'M','W','D'}
        date_format='mmmyy';
end
%       10             'yyyy'                   2000
%       11             'yy'                     00
%       12             'mmmyy'                  Mar00
%       17             'QQ-YY'                  Q1-96
%       28             'mmmyyyy'                Mar2000

pp.xtick_labels=datestr(pp.tickLocs,date_format);
pp.xlim=[pp.xdatenums(1),pp.xdatenums(end)];