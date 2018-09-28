classdef rise_file < rise_report.generic_report
    
    properties
        % type face applies to Captions only...
        typeface ={'\large','\bfseries'} % \textbf,\textit,\underline,\emph,\large, etc. or a combination
        
        filename = ''
        
        latexAlias = false
        
        centering=false
        
        lines = @all
        
        lineNumbers=true
        
        footnote=''
        
        syntax = true
        
        default={ ...
            'latexalias',false,@islogicalscalar,false, ...
            'linenumbers',true,@islogicalscalar,true, ...
            'lines',@all,@(x) isequal(x,@all) || isnumeric(x),true, ...
            'syntax',true,@(x)islogical && isscalar,true
            }
        
    end
    
    properties(Dependent)
        
        batch
        
    end
    
    
    methods
        
        function obj = rise_file(varargin)% title,filename,varargin
            
            obj=rise_report.feed_properties(mfilename,obj,varargin{:});
            
        end
        
        function b=get.batch(obj)
            
            opts=struct('lines',obj.lines,'syntax',obj.syntax,...
                'linenumbers',obj.lineNumbers,'latexalias',obj.latexAlias);
            
            b= print_model_engine(obj.filename,opts);
            
            b = reformat_text(obj,b);
            
            b=regexp(b,'\n','split');
            
        end
        
    end
    
end


function C = print_model_engine(filename,options)
% printmodelfile  [Not a public function] LaTeXify and syntax highlight model file.
%
rise_keywords=load_keywords();

ckwd=parser.cell2matize(rise_keywords);

C = '';

if isempty(filename)
    
    return
    
end

opts={
    'lines',@all
    'syntax',true
    'linenumbers',true
    'latexalias',true
    };

if nargin<2
    
    options=[];
    
end

if isempty(options)
    
    options=cell2struct(opts(:,2),opts(:,1),1);
    
end

br = sprintf('\n'); %#ok<SPRINTFN>

% Read the text file into a cellstr with EOLs removed.
txt = parser.read_file(filename,false);

if isequal(options.lines,@all)
    
    numberOfLines = txt{end,3};
    
    options.lines = 1 : numberOfLines;
    
else
    
    txt = txt(options.lines);
    
    numberOfLines = size(txt,1);
    
end

% Choose escape character.
escList = '`@?$#~&":|!^[]{}<>';

esc = pick_one_escape_character(txt(:,1),escList);

if isempty(esc)
    
    utils.error('report', ...
        ['Cannot print the model file. ', ...
        'Make sure at least on of these characters completely disappears ', ...
        'from the model file: %s.'], ...
        escList);
    
end

verbEsc = ['\verb',esc];

nDigit = ceil(log10(max(options.lines)));

C = [C,'\definecolor{mylabel}{rgb}{0.55,0,0.35}',br];

C = [C,'\definecolor{myparam}{rgb}{0.90,0,0}',br];

C = [C,'\definecolor{mykeyword}{rgb}{0,0,0.75}',br];

C = [C,'\definecolor{mycomment}{rgb}{0,0.50,0}',br];

quote_open=false;

close_quote='"__';

open_quote='__"';

for iline = 1 : numberOfLines
    
    c = fileOneLine(txt{iline,1});
    
    C = [C,c,' \\',br]; %#ok<AGROW>
    
end

C = strrep(C,['\verb',esc,esc],'');

C=strrep(C,open_quote,[esc,'{\color{myparam}',verbEsc,'"']);

C=strrep(C,close_quote,['"',esc,'}',verbEsc]);

C=strrep(C,'?',[esc,'{\color{myparam}\',verbEsc,'?',esc,'}']);

kwd3=parser.cell2matize({'#','!','?'});

attavila=[esc,'{\\color{myparam}\',verbEsc,'$1',esc,'}\',verbEsc,'$2'];

C=regexprep(C,[kwd3,'(\s*\w+\s*=)'],attavila);

    function C = fileOneLine(C)
        
        lineComment = '';
        
        pos = strfind(C,'%');
        
        if ~isempty(pos)
            
            pos = pos(1);
            
            lineComment = C(pos:end);
            
            C = C(1:pos-1);
            
        end
        
        if options.syntax
            % Keywords.
            
            keywordsFunc = @doKeywords; %#ok<NASGU>
            
            C = regexprep(C, ['\<',ckwd,'\>'], ...
                '${keywordsFunc($0)}');
            
            % Line comments.
            if ~isempty(lineComment)
                
                lineComment = [ ...
                    esc, ...
                    '{\color{mycomment}', ...
                    verbEsc,lineComment,esc, ...
                    '}', ...
                    verbEsc];
                
            end
            
            % double quotes
            tmpC='';
            
            while true
                
                bingo=find(C=='"',1,'first');
                
                if isempty(bingo)
                    
                    break
                
                end
                
                if quote_open
                    % time to close
                    str=close_quote;
                    
                else
                    % time to open
                    str=open_quote;%'__"';
                    
                end
                
                tmpC=[tmpC,C(1:bingo-1),str]; %#ok<AGROW>
                
                C=C(bingo+1:end);
                
                quote_open=~quote_open;
                
            end
            
            tmpC=[tmpC,C];
            
            C=tmpC;
            
        end
        
        if options.linenumbers
            
            C = [ ...
                sprintf('%*g: ',nDigit,options.lines(iline)), ...
                C];
            
        end
        
        C = [verbEsc,C,lineComment,esc];
        
        function C = doKeywords(C)
            
            color = 'mykeyword';
            
            C = [ ...
                '{\color{',color,'}', ...
                verbEsc,C,esc, ...
                '}', ...
                ];
            C = [esc,C,verbEsc];
        end
        
    end 

end

function Esc = pick_one_escape_character(File,EscList)

File = [File{:}];

Esc = '';

for i = 1 : length(EscList)
    
    if ~parser.mycontains(File,EscList(i))
        
        Esc = EscList(i);
        
        break
        
    end
    
end

end 

function kwd=load_keywords()

rise_keywords=parser.initialize_blocks();

rise_keywords=rise_keywords(:,2);

kwd=cell(1,1000);

iter=0;

for ii=1:numel(rise_keywords)
    
    kwdi=rise_keywords{ii};
    
    if ischar(kwdi)
        
        kwdi={kwdi};
        
    end
    
    n=numel(kwdi);
    
    kwd(iter+(1:n))=kwdi(:).';
    
    iter=iter+n;
    
end

kwd=kwd(1:iter);

kwd2=strcat('@#\s*',{'if','else','end','switch','case','for','while','stst','sstate'});

kwd=[kwd,kwd2];

end


