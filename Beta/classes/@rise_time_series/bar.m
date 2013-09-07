function [plot_handle,axis_handle]=bar(this,varargin)
datta=double(this);
if size(datta,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
date_numbers=vertcat(this.TimeInfo.date_number);
tmp=bar(date_numbers,datta,varargin{:});
NumTicks = min(5,numel(date_numbers));
tick_locs=round(linspace(1,100,NumTicks)/100*numel(date_numbers));
tick_locs(tick_locs==0)=1;
tick_locs=date_numbers(tick_locs);
switch this.frequency
    case ''
        date_format='yy';
    case {'H','Q'}
        date_format='QQ-YY';
    case {'M','W','D'}
        date_format='mmmyy';
end

tick_labels=datestr(tick_locs,date_format);
set(gca,'XLim',[date_numbers(1),date_numbers(end)],...
    'XTick',tick_locs,...
    'XTickLabel',tick_labels)
%       10             'yyyy'                   2000
%       11             'yy'                     00
%       12             'mmmyy'                  Mar00
%       17             'QQ-YY'                  Q1-96
%       28             'mmmyyyy'                Mar2000


if nargout>0
    plot_handle=tmp;
    if nargout>1
        axis_handle=get(gca);
    end
end
end

% function H=bar(this,varargin)
% date_numbers=vertcat(this.TimeInfo.date_number);
% datta=double(this);
% if size(datta,3)>1
%     error([mfilename,':: this operation is only defined for databases with one page'])
% end
% tmp=bar(date_numbers,datta,varargin{:});
% NumTicks = 8;
% tick_locs=round(linspace(1,100,NumTicks)/100*numel(date_numbers));
% tick_locs(tick_locs==0)=1;
% tick_locs=date_numbers(tick_locs);
% tick_labels=datestr(tick_locs,11);
% set(gca,'XLim',[date_numbers(1),date_numbers(end)],...
%     'XTick',tick_locs,...
%     'XTickLabel',tick_labels)
% legend(this.varnames)
% if nargout>0
%     H=tmp;
% end
% end
