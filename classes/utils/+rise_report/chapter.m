classdef chapter < rise_report.generic_report
    properties
        title=''
        tableOfContents=false
        tableOfContentsTitle=''
        numbering=true % preface is a chapter and is not numbered
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=chapter(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
        end
        function b = get.batch(obj)
            b=rise_report.title_item(mfilename,obj.title,obj.numbering,...
                obj.tableOfContentsTitle);
        end
    end
end
