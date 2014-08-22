function [h,pp]=plot_real_time(rts,pp,nticks)
if nargin<3
    nticks=[];
end
dd=double(rts);
ymax=max(max(dd));
ymin=min(min(dd));
[nr,nc]=size(dd);
if nargin<2
    pp=[];
end
if isempty(pp)
    xtime=rts.date_numbers(1)-1:rts.date_numbers(end)+nc;
    pp=plot_specs(xtime,nticks);
end

map=getappdata(0,'rise_default_plot_colors');
plot(pp.xdatenums(1:nr)',dd(:,1),'linewidth',2)%,'color',[0,0,0]
hold on
iter=0;
for tt=1:size(dd,1)
    iter=iter+1;
    plot(pp.xdatenums(tt+(0:nc-1))',dd(tt,:)','color',map{iter})%,'linestyle','--'
    if iter==numel(map)
        iter=0;
    end
end

set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels,...
    'ylim',[ymin,ymax])
hold off
grid on

h=gca;

end
