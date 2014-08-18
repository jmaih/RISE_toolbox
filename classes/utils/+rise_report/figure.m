classdef figure < rise_report.generic_report
    properties
        title=''
        % string or function handle
        name ='' 
        % scale graphic to the specified width
        width
        % scale graphic to the specified height
        height
        % rotate graphic counterclockwise
        angle=0
        % scale graphic
        scale=0.85
        % path to the graph
        graphicspath = ''
        numbering=true
    end
    properties(Dependent)
        batch
    end
    properties(Hidden)
        figure_number
    end
    methods
        function obj=figure(varargin)
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            own_props={
                'name',@(x)ischar(x)||ishandle(x)
                'title',@(x)ischar(x)
                'width',@(x)isnumeric(x) && x<=1 && x>0
                'height',@(x)isnumeric(x) && x<=1 && x>0
                'angle',@(x)isnumeric(x) && x<=360 && x>=-360
                'scale',@(x)isnumeric(x) && x<=1 && x>0
                'numbering',@(x)islogical(x)
                };
            nop=size(own_props,1);
            for op=1:nop
                if ~isempty(obj.(own_props{op,1}))
                    assert(own_props{op,2}(obj.(own_props{op,1})))
                end
            end
        end
        function b = get.batch(obj)
            b={};
            if ~isempty(obj.name)
                titel=obj.title;
                if ~isempty(titel)
                    if obj.numbering
                        titel=sprintf('Figure \\# %0.0f: %s',obj.figure_number,titel);
                    end
                    reprocess=@(x)rise_report.generic_report.reprocess(x);
                    % best_title=@(x)rise_report.generic_report.best_title(x);
                    %                     b=[b
                    %                         '\centering'
                    %                         best_title(obj.title)];
                end
                tmpname=obj.name;
                angle_=obj.angle;
                if ishandle(tmpname)
                    new_handle=tmpname;
                    tmpname=sprintf('%s/figure_%0.0f',...
                        obj.graphicspath,obj.figure_number);% tmpname=tempname(obj.graphicspath);
                    angle_=utils.plot.saveaspdf(new_handle,tmpname);
                end
                tmpname=strrep(tmpname,'\','/');%tmpname=reprocess(strrep(tmpname,'\','/'));
                tmpname=parser.remove_file_extension(tmpname);%,'.pdf'
                % quotes and spaces may still not work if the extension .pdf is
                % added
                if any(isspace(tmpname))
                    tmpname=['"',tmpname,'"'];
                end
                myfigname=reprocess(titel);
                attributes={};
                if ~isempty(obj.width)
                    attributes=[attributes,sprintf('%0.2f\textwidth',obj.width)];
                end
                if ~isempty(obj.height)
                    attributes=[attributes,sprintf('%0.2f\textheight',obj.height)];
                end
                if ~isempty(angle_)
                    attributes=[attributes,sprintf('angle=%0.2f',angle_)];
                end
                if ~isempty(obj.scale)
                    attributes=[attributes,sprintf('scale=%0.2f',obj.scale)];
                end
                attributes=cell2mat(strcat(attributes,','));
                b=[b
                    '\begin{tabular}[H]{@{\hspace*{-3pt}}c@{}}'
                    ['\multicolumn{1}{c}{\large\bfseries ',myfigname,'}\\']
                    ['\centerline{\includegraphics[',attributes(1:end-1),...
                    ']{',tmpname,'}}']
                    '\end{tabular}'
                    ];
            end
        end
    end
end