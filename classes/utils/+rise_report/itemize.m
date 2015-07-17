classdef itemize < rise_report.generic_report
    % itemize report itemize object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.itemize/addlistener)
    % - [best_title](rise_report.itemize/best_title)
    % - [delete](rise_report.itemize/delete)
    % - [eq](rise_report.itemize/eq)
    % - [findobj](rise_report.itemize/findobj)
    % - [findprop](rise_report.itemize/findprop)
    % - [ge](rise_report.itemize/ge)
    % - [gt](rise_report.itemize/gt)
    % - [isvalid](rise_report.itemize/isvalid)
    % - [itemize](rise_report.itemize/itemize)
    % - [le](rise_report.itemize/le)
    % - [lt](rise_report.itemize/lt)
    % - [ne](rise_report.itemize/ne)
    % - [notify](rise_report.itemize/notify)
    % - [reprocess](rise_report.itemize/reprocess)
    % - [write](rise_report.itemize/write)
    %
    % properties
    % -----------
    %
    % - [items] -
    % - [markers] -
    % - [batch] -
    % - [id] -
    properties
        items={} % cellstring
        % cellstrings to customize the things of the items e.g.
        % \item[Stupid] \item[Smart] \item[-]
        markers={}
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=itemize(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            if ~isempty(obj.markers) && numel(obj.markers)~=numel(obj.items)
                if numel(obj.markers)==1
                    obj.markers=repmat(obj.markers,size(obj.items));
                else
                    error('number of markers different from number of items')
                end
            end
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.items)
                marker='';
                marker_style=~isempty(obj.markers);
                b={['\begin{',mfilename,'}']};
                for it=1:numel(obj.items)
                    if marker_style
                        marker=['[',obj.markers{it},']'];
                    end
                    b=[b
                        ['\item',marker,' ',obj.items{it}]
                        ]; %#ok<AGROW>
                end
                b=[b
                    ['\end{',mfilename,'}']
                    ];
            end
        end
    end
end

