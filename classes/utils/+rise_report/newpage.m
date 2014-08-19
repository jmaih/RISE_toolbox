classdef newpage < rise_report.generic_report
    properties(Constant)
        % plain prints the page numbers on the bottom of the page, in the
        % middle of the footer. This is the default page style.
        % headings prints the current chapter heading and the page number
        % in the header on each page, while the footer remains empty. (This
        % is the style used in this document)
        % empty sets both the header and the footer to be empty.
        thispagestyle=''
    end
    properties(Dependent)
        batch
    end
    methods
        function obj=newpage(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            if nargin
                if ~any(strcmp(varargin{2},rise_report.page_styles()))
                    error([varargin{2},' is not a valid page style'])
                end
            end
        end
        function b = get.batch(obj)
            b={' ';'\newpage'};
            if ~isempty(obj.thispagestyle)
                b=[
                    b
                    ['\thispagestyle{',obj.thispagestyle,'}']
                    ];
            end
        end
    end
end
