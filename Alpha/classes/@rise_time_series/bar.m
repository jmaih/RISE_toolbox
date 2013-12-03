function [plot_handle,axis_handle]=bar(this,varargin)
datta=double(this);
if size(datta,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end

nticks=[];
n=length(varargin);
extract=[];
for icol=1:n
    if strcmpi(varargin{icol},'nticks')
        nticks=varargin{icol+1};
        extract=[icol,icol+1];
        break
    end
end
if ~isempty(extract)
    varargin(extract)=[];
end

pp=plot_specs(this.date_number,nticks);

tmp=bar(pp.xdatenums,datta,varargin{:});

set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels) %...

grid on

if nargout>0
    plot_handle=tmp;
    if nargout>1
        axis_handle=get(gca);
    end
end
end