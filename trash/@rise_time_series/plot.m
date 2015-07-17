function [plot_handle,axis_handle]=plot(this,varargin)
datta=double(this);
if size(datta,2)>1 && size(datta,3)>1
    error([mfilename,':: cannot handle many variables and many pages simultaneously'])
else
    datta=squeeze(datta);
end

[nticks,varargs]=extract_number_of_ticks(varargin{:});

pp=plot_specs(this.date_number,nticks);
% %--------------------
% [stamp,unstamp]=time_frequency_stamp();
% s=this.date_number;
% freq=unstamp(s(1)-floor(s(1)));
% year=floor(s/freq);
% period=round(s-year*freq+1-stamp(freq));
% % convert period to months
% switch freq
%     case 1
% dateFormat =10; %<--'yyyy';
%     case {2,4}
% dateFormat =17; %<--'yyyy';
%     case 12
% dateFormat =12; %<--'yyyy';
%     otherwise
% end
% if freq~=1
%     period=period*12/freq;
% end
% dd=datenum(year,period,1);
% %--------------------
% tmp=plot(dd,datta,varargin{:});
% datetick('x',dateFormat);
tmp=plot(pp.xdatenums,datta,varargs{:});
set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %,'XTickMode','auto'...
%     'ylim',[ymin,ymax])

grid on

if nargout>0
    plot_handle=tmp;
    if nargout>1
        axis_handle=get(gca);
    end
end
end

function [plot_handle]=line(this,varargin)
date_numbers=this.date_number;
datta=double(this);
if size(datta,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
NumTicks = 8;
tick_locs=round(linspace(1,100,NumTicks)/100*numel(date_numbers));
tick_locs(tick_locs==0)=1;
tick_locs=date_numbers(tick_locs);
tick_labels=datestr(tick_locs,11);
if this.NumberOfVariables==2
    line(date_numbers,datta(:,1),'Color','r');
    ax1 = gca;
    %                 set(ax1,'XColor','r','YColor','r')
    ax2 = axes('Position',get(ax1,'Position'),...
        'YAxisLocation','right',...
        'Color','none');
    %                     'XAxisLocation','top',...
    %                     'XColor','k','YColor','k'
    line(date_numbers,datta(:,2),'Color','k','Parent',ax2);
    set([ax1,ax2],'XLim',[date_numbers(1),date_numbers(end)],...
        'XTick',tick_locs,...
        'XTickLabel',tick_labels)
else
    line(date_numbers,datta,varargin{:});
    set(gca,'XLim',[date_numbers(1),date_numbers(end)],...
        'XTick',tick_locs,...
        'XTickLabel',tick_labels)
end
if nargout>0
    plot_handle=gca;
end
end

