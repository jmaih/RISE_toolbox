function pp=plot_specs(date_numbers,nticks,date_format)
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

default_nticks=8;
if nargin<3
    date_format=[];
    if nargin<2
        nticks=[];
        if nargin<1
            pp=default_nticks;
            return
        end
    end
end
if isempty(nticks)
    nticks=default_nticks;
end
% see also plot_real_time

pp=struct();

pp.xdatenums=date_numbers;

pp.nticks=min(nticks,numel(pp.xdatenums));

tickLocs=utils.plot.locate_ticks(numel(pp.xdatenums),pp.nticks);

pp.tickLocs=pp.xdatenums(tickLocs);

[pp.xtick_labels,frequency0]=serial2date(pp.tickLocs);

if ~isempty(date_format)
    [y,p0]=date2year_period(pp.tickLocs);
    freq0=frequency2num(frequency0);
    freq1=12;
    head=false;
    p1=period2period(p0(:),freq0,freq1,head);
    dv=[y(:),p1(:)];
    dv=[dv,ones(size(dv,1),1),zeros(size(dv,1),3)];
    pp.xtick_labels=datestr(dv,date_format);
end

pp.xlim=[pp.xdatenums(1),pp.xdatenums(end)];
% set(gca,'xtickLabel',{'alpha','bravo','charlie'})
end