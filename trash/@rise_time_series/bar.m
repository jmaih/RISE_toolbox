function [plot_handle,axis_handle]=bar(this,varargin)
datta=double(this);
if size(datta,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end

[nticks,varargs]=extract_number_of_ticks(varargin{:});

pp=plot_specs(this.date_number,nticks);

tmp=bar(pp.xdatenums,datta,varargs{:});

set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %...

grid on

if nargout>0
    plot_handle=tmp;
    if nargout>1
        axis_handle=get(gca);
    end
end
end