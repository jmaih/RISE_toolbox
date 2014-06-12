function plot_packages(db,r0,c0,packages,start,finish,figtitle,horizline,func_plot,tight)

if nargin<10
    tight=[];
if nargin<9
    func_plot=[];
    if nargin<8
        horizline=[];
    end
end
end
if isempty(tight)
    tight=false;
end
if isempty(r0)
    r0=3;
end
if isempty(c0)
    c0=3;
end

if isa(db,'ts')
    db=pages2struct(db);
end
if isempty(start)
    vnames=fieldnames(db);
    start=db.(vnames{1}).start;
end
if isempty(finish)
    vnames=fieldnames(db);
    finish=db.(vnames{1}).finish;
end

[zz,zz_name]=get_zero_time_series();
zline=~isempty(zz);

nvar=numel(packages);
nstar=r0*c0;
nfigs=ceil(nvar/nstar);
str='';
thisFuncPlot=func_plot;
for fig=1:nfigs
    if nfigs>1
        str=['(',sprintf('%0.0f',fig),')'];
    end
    [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,nvar,r0,c0);
    figure('name',[figtitle,str])
    for iplot=1:min(r*c,Remains)
        iobs=(fig-1)*nstar+iplot;
        vname=packages{iobs};
        char_flag=ischar(vname);
        subplot(r,c,iplot)
        if isempty(func_plot)
            thisFuncPlot=@plot;
        end
        if char_flag
            tmp=[db.(vname),zz];
        else
            if numel(vname)==2 && ~zline && isempty(func_plot)
                thisFuncPlot=@plotyy;
            end
            tmp=db.(vname{1});
            for iname=2:numel(vname)
                tmp=[tmp,db.(vname{iname})]; %#ok<*AGROW>
            end
            tmp=[tmp,zz];
            tmp.varnames=[vname,zz_name];
        end
        plot_window(tmp,start,finish,thisFuncPlot,'linewidth',2);
        if char_flag
            title(vname,'interp','none')
        else
            legend(vname)
        end
        if zline
            zzdata=double(zz);
            set(findobj(gca,'YData',zzdata),'color',[1,0,0])
        end
        if tight
            axis tight;
        else
            tmp=double(tmp);
            lowest=min(tmp(:));
            highest=max(tmp(:));
            ymin=lowest-abs(lowest/10);
            ymax=highest+abs(highest/10);
            ylim([ymin,ymax])
        end
    end
end


    function [zz,zz_name]=get_zero_time_series()
        zz=ts.empty(1,0);
        zz_name='';
        if ~isempty(horizline)
            zz_name='horizline';
            start_serial=date2serial(start);
            finish_serial=date2serial(finish);
            nobs=numel(start_serial:finish_serial);
            zz=ts(start,horizline*ones(nobs,1),zz_name);
        end
    end
end

%{
color codes
RGB Value, Short Name Long Name
[1 1 0] y yellow
[1 0 1] m magenta
[0 1 1] c cyan
[1 0 0] r red
[0 1 0] g green
[0 0 1] b blue
[1 1 1] w white
[0 0 0] k black
%}