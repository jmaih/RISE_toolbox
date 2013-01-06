function [h,pp]=plot_real_time(rts,pp)
dd=double(rts);
ymax=max(max(dd));
ymin=min(min(dd));
[nr,nc]=size(dd);
if nargin<2
    pp=[];
end
if isempty(pp)
    xtime=rts.TimeInfo(1)-1:rts.TimeInfo(end)+nc;
    nticks=8;
    pp=plot_specs(xtime,nticks);
end

map=getappdata(0,'rise_default_plot_colors');
plot(pp.xdatenums(2:nr+1)',dd(:,1),'linewidth',2)%,'color',[0,0,0]
hold on
iter=0;
for tt=1:size(dd,1)
    iter=iter+1;
    plot(pp.xdatenums(tt+(0:nc-1))',dd(tt,:)','color',map{iter})%,'linestyle','--'
    if iter==numel(map)
        iter=0;
    end
end
% add the steady state
% vloc=locate_variables(obs_names{iobs},{lwz.varendo.name});
% plot(pp.xlim,steady_state*ones(1,2),'color',[1,0,0])
set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels,...
    'ylim',[ymin,ymax])
hold off
grid on

h=gca;

end
