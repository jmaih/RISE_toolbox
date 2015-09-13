function [tex_name,legend_]=one_fanchart(f,vname,mycol,histdb,obsnames)
data=f.(vname);
plot_fanchart(data,mycol);
if any(strcmp(vname,obsnames))
    S=struct('type','()','subs',{{vname}});
    obj=subsref(histdb,S);
    hold on
    plot('2005:2010',obj,'linewidth',2)
    hold off
end
legend_=[];
tex_name=vname;
end