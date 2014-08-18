classdef verbatim < rise_report.text
    methods
        function obj=verbatim(varargin)
            obj=obj@rise_report.text(varargin{:});
        end
    end
    methods(Access = private)
        function b=batch_implementation(~,b)
            b=['\begin{verbatim}';b(:);'\end{verbatim}'];
        end
    end
end
