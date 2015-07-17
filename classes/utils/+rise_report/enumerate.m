classdef enumerate < rise_report.generic_report
    % enumerate report enumerate object
    %
    % methods
    % --------
    %
    % - [addlistener](rise_report.enumerate/addlistener)
    % - [best_title](rise_report.enumerate/best_title)
    % - [delete](rise_report.enumerate/delete)
    % - [enumerate](rise_report.enumerate/enumerate)
    % - [eq](rise_report.enumerate/eq)
    % - [findobj](rise_report.enumerate/findobj)
    % - [findprop](rise_report.enumerate/findprop)
    % - [ge](rise_report.enumerate/ge)
    % - [gt](rise_report.enumerate/gt)
    % - [isvalid](rise_report.enumerate/isvalid)
    % - [le](rise_report.enumerate/le)
    % - [lt](rise_report.enumerate/lt)
    % - [ne](rise_report.enumerate/ne)
    % - [notify](rise_report.enumerate/notify)
    % - [reprocess](rise_report.enumerate/reprocess)
    % - [write](rise_report.enumerate/write)
    %
    % properties
    % -----------
    %
    % - [items] -
    % - [markers] -
    % - [batch] -
    % - [id] -
    properties % cellstring
        items={}
        % cellstrings to customize the things of the items e.g.
        % \item[Stupid] \item[Smart] \item[-]
        markers={}
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=enumerate(varargin)
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

