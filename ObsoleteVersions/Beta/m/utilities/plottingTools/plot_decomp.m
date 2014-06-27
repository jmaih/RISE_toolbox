function h=plot_decomp(obj)

uall=double(obj);
upos=zeros(size(uall));
uneg=zeros(size(uall));
ipos=uall>=0;
ineg=uall<0;
upos(ipos)=uall(ipos);
uneg(ineg)=uall(ineg);
this_pos=rise_time_series(obj.TimeInfo,upos,obj.varnames);
this_neg=rise_time_series(obj.TimeInfo,uneg,obj.varnames);
this_all=rise_time_series(obj.TimeInfo,sum(uall,2));
bar(this_pos,'stack') %  area(this_pos)
hold on
bar(this_neg,'stack') % area(this_neg)
hold on
plot(this_all,'k-','linewidth',2)
axis tight;
hold off
if nargout>0
    h=gca;
end
