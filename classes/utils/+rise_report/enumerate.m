classdef enumerate < rise_report.generic_report
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

