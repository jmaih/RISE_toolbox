function pp=plot_specs(date_number,nticks)
if nargin<2
    nticks=[];
end
if isempty(nticks)
    nticks=8;
end
% see also plot_real_time

pp=struct();
pp.xdatenums=date_number;
pp.nticks=min(nticks,numel(pp.xdatenums));
tickLocs=ceil(linspace(1,numel(pp.xdatenums),pp.nticks));
tickLocs=unique(tickLocs); % just in case I accidentally add the same element twice from the lines above
% discard the even locations
pp.tickLocs=pp.xdatenums(tickLocs);

% [~,frequency]=serial2date(date_number(1));
% switch frequency
%     case ''
%         date_format='yy';
%     case {'H','Q'}
%         date_format='QQ-YY';
%     case {'M','W','D'}
%         date_format='mmmyy';
% end
%       10             'yyyy'                   2000
%       11             'yy'                     00
%       12             'mmmyy'                  Mar00
%       17             'QQ-YY'                  Q1-96
%       28             'mmmyyyy'                Mar2000

pp.xtick_labels=serial2date(pp.tickLocs);% datestr(pp.tickLocs,date_format);
pp.xlim=[pp.xdatenums(1),pp.xdatenums(end)];
% set(gca,'xtickLabel',{'alpha','bravo','charlie'})