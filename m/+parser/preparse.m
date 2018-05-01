function [output,has_macro]=preparse(FileName,parsing_definitions)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% {'movave','movavg','movgeom','movprod','movsum','diff','difflog'}

% @#include(filename)
% @#include"filename"
% @#import(filename)
% @#import"filename"

% @#for co={ca,fr,de,be}
%     Y_@{co}=rho_@{co}*Y_@{co}{-1}+sig_@{co}*E_@{co};
%     @#end
% 
% @#for co=[1,4,6,8]
%     Y_@{co}=rho_@{co}*Y_@{co}{-1}+sig_@{co}*E_@{co};
%     @#end
% 
% @#for co=1:8
%     Y_@{co}=rho_@{co}*Y_@{co}{-1}+sig_@{co}*E_@{co};
%     @#end
% 
% @#for co=1:someNumber (where someNumber is defined outside...)
%     Y_@{co}=rho_@{co}*Y_@{co}{-1}+sig_@{co}*E_@{co};
%     @#end
% 
% @#for co=1:someNumber (where someNumber is defined outside...)
%     Y_@{co}=rho_@{co}*Y_@{@{co}+2}{-1}; % evaluate @{co} then @{@{co}+2}
%     @#end
% 
% @#for co=1:8
%     @#for coo={ca,fr,de,be}
%         Y_@{co}_@{coo}=rho_@{co}_@{coo}*Y_@{co}_@{coo}{-1}+sig_@{co}_@{coo}*E_@{co}_@{coo};
%         @#end
%     @#end

% @# if strcmp(capRec,'firstCase')
%     @# elseif strcmp(capRec,'anotherCase')
%     @# end

% @# switch capRec
%     @# case 'firstCase'
%     @# case 'secondCase'
%     @# otherwise
%     @# end

% int2str(x)=sprintf('%.0f',x)
if nargin<2
    
    parsing_definitions=struct();
    
    if nargin<1
        
        output=struct('rise_flags',struct());
        
        return
        
    end
    
end

if ~isstruct(parsing_definitions)
    
    if ~iscell(parsing_definitions)||(iscell(parsing_definitions) && size(parsing_definitions,2)~=2)
        
        error('parsing definitions must be a struct or a cell array with two columns')
        
    end
    
    tmp=struct();
    
    for irow=1:size(parsing_definitions,1)
        
        tmp.(parsing_definitions{irow,1})=parsing_definitions{irow,2};
        
    end
    
    parsing_definitions=tmp;
    
end

% first sweep: read the file, and remove comments and stamp the lines
RawFile=parser.read_file(FileName);

% second sweep: process the if else end for and include comments
% include should also go through preparsing
[output,has_macro]=preparse_expand(RawFile,parsing_definitions);

end

%--------------------------------------------------------------------------

function [output,has_macro]=preparse_expand(rawfile,definitions,has_macro,output)
if nargin<4
    
    output=[];
    
    if nargin<3
        
        has_macro=false;
        
    end
    
end

if isempty(rawfile)
    
    return
    
end

remove_space=@(x)x(~isspace(x));

if isempty(output),output=cell(0,3);end

if isempty(has_macro),has_macro=false;end

has_macro=false;

process_include()

[rawfile,has_macro]=process_flow_controls(rawfile,definitions,has_macro);

% evaluate remaining @{function} if possible
%--------------------------------------------
myevaluate=@reevaluate; %#ok<NASGU>

rawfile(:,1)=regexprep(rawfile(:,1),'(@\{[^\}]+\})','${myevaluate($1)}');

% % process keywords last
% rawfile(:,1)=process_keywords(rawfile(:,1));

output=[output;rawfile];

    function process_include()
        
        arobase=include_arobase(rawfile(:,1));
        
        while ~isempty(arobase)
            
            if ~has_macro
                
                has_macro=true;
                
            end
            
            process_include_engine()
            
        end
        
        function process_include_engine
            
            tmp=rawfile(1:0,:);
            
            locs=[arobase.pos];
            
            offset=0;
            
            for iloc=1:numel(locs)
                
                rline_=next_line(locs(iloc));
                
                tmp=[tmp
                    rawfile(offset+1:locs(iloc)-1,:);
                    ]; %#ok<AGROW>
                
                tmp=parser.append_file(tmp,rline_,definitions);
                
                offset=locs(iloc);
                
            end
            
            tmp=[tmp
                rawfile(offset+1:end,:);
                ];
            
            rawfile=tmp;
            
            arobase=include_arobase(rawfile(:,1));
            
        end
        
    end

    function rline=next_line(iter)
        
        rline.rawline = rawfile{iter,1};
        
        rline.FileName= rawfile{iter,2};
        
        rline.row_number_string= int2str(rawfile{iter,3});
        
        [startIndex,endIndex] = regexp(rline.rawline,['@\s*(?#0 or more spaces)',...
            '#\s*(include|if|for|else(if)?|end(if|for)?',...
            '(?#optionally followed by if or for))']);
        
        if numel(startIndex)>1
            
            error('')
            
        elseif numel(startIndex)==1
            
            rline.first_tok=remove_space(rline.rawline(startIndex:endIndex));
            
            rline.rest=rline.rawline(endIndex+1:end);
            
        else
            
            rline.first_tok='';
            
            rline.rest=rline.rawline;
            
        end
        
    end

    function string=reevaluate(string)
        
        left=find(string=='{');
        
        right=find(string=='}');
        
        middle=string(left+1:right-1);
        
        try
            
            middle=eval_action(definitions,middle);
            
            if ~isnan(middle)
                
                string=sprintf('%0.15g',middle);
                
            end
            
        catch me
            
            warning(me.message)
            
        end
        
    end

end

%--------------------------------------------------------------------------

function [rawfile,has_macro]=process_flow_controls(rawfile,definitions,has_macro)

arobase=[];

update_arobase();

while ~isempty(arobase)
    
    if ~has_macro
        
        has_macro=true;
        
    end
    
    process_flow_control_engine()
    
end

    function process_flow_control_engine()
        
        locs=[arobase.pos];
        
        switch arobase(1).type
            case 'for'
                do_for()
            case 'if'
                do_if()
            case 'switch'
                do_switch()
        end
        
        update_arobase()
        
        function do_for()
            
            top=rawfile(1:locs(1)-1,:);
            
            bottom=rawfile(locs(end)+1:end,:);
            
            middle=rawfile(locs(1)+1:locs(end)-1,:);
            
            info=regexp(arobase(1).action,'(?<stud>\w+)\s*(=|in)\s*(?<set>.+)','names');
            
            update_set()
            
            for ii=1:numel(info.set)
                
                tmp=middle;
                
                tmp(:,1)=regexprep(tmp(:,1),['@{',info.stud,'}'],info.set{ii});
                
                top=[top;tmp]; %#ok<AGROW>
                
            end
            
            rawfile=[top;bottom];
            
            function update_set()
                
                sset=info.set;
                
                if any(sset=='{')||any(sset=='[')
                    
                    sset=regexp(sset,',','split');
                    
                    sset=regexprep(sset,'(\{|\}|\[|\])','');
                    
                elseif any(sset==':')
                    
                    sset=regexp(sset,':','split');
                    
                    if numel(sset)==2
                        
                        sset=[sset(1),'1',sset(2)];
                        
                    end
                    
                    for iii=1:numel(sset)
                        
                        if ~all(isstrprop(sset{iii},'digit'))
                            
                            sset{iii}=eval_action(definitions,sset{iii});
                            
                            sset{iii}=sprintf('%0.0g',sset{iii});
                            
                        end
                        
                    end
                    
                    sset=str2double(sset{1}):str2double(sset{2}):str2double(sset{3});
                    
                    sset=num2cell(sset);
                    
                    sset=cellfun(@(x)num2str(x),sset,'uniformOutput',false);
                    
                else % evaluate the set
                    
                    sset=eval_action(definitions,sset);
                    
                    if isa(sset,'double')
                        
                        sset=num2cell(sset);
                        
                    end
                    
                    % the set is not necessarily a cellstr and so we have
                    % to apply some checking/correction
                    
                    if ~iscellstr(sset)
                        
                        sset=cellfun(@(x)num2str(x),sset,'uniformOutput',false);
                        
                    end
                    
                end
                
                info.set=sset;
                
            end
            
        end
        
        function do_if()
            
            start=1;
            
            ff=@(x)eval_action(definitions,x);
            
            rmQuotes=@(x)x;
            
            do_if_or_switch(start,ff,rmQuotes)
            
        end
    
        function do_switch()
            
            test=regexprep(rawfile(locs(1)+1:locs(2)-1,1),'%.+','');
            
            if ~all(cellfun(@isempty,test,'uniformOutput',true))
                
                error('Non comment found between switch and case')
                
            end
            
            stud=arobase(1).action;
            
            start=2;
            
            ff=@(x)strcmp(definitions.(stud),x);
            
            rmQuotes=@(x)strrep(x,'''','');
            
            do_if_or_switch(start,ff,rmQuotes)
        
        end
        
        function do_if_or_switch(start,ff,rmQuotes)
            
            for icase=start:numel(arobase)-1
                
                action=rmQuotes(arobase(icase).action);
                
                action(isspace(action))=[];
                
                if isempty(action)
                    % otherwise or else
                    found=true;
                    
                else
                    
                    found=ff(action);
                    
                end
                
                if found
                    
                    break
                    
                end
                
            end
            
            if found
                
                discard=[locs(1):locs(icase),locs(icase+1):locs(end)];
                
            else
                
                discard=locs(1):locs(end);
                
            end
            
            rawfile(discard,:)=[];
            
        end
        
    end

    function update_arobase()
        
       [arobase]=control_flow_arobase(rawfile(:,1));
        
    end

end

%--------------------------------------------------------------------------

function flag___=eval_action(def____,action)

fields___=fieldnames(def____);

for ifield___=1:numel(fields___)
    
    name______=fields___{ifield___};
    
    eval([name______,'=def____.(name______);'])
    
end

flag___=eval(action);

end

%--------------------------------------------------------------------------

function [arobase]=include_arobase(txt)

% @#include(filename)
% @#include"filename"
% @#import(filename)
% @#import"filename"

arobase=regexp(txt,'@#\s*(include|import)\s*("|\()\s*(?<filename>[^"\)]+)("|\))(?<pos>)','names');

for ib=1:numel(arobase)
    
    if ~isempty(arobase{ib})
        
        arobase{ib}.pos=ib;
        
    end
    
end

arobase=[arobase{:}];

end

%--------------------------------------------------------------------------

function [arobase]=control_flow_arobase(txt)

triggers={'switch','if','for'};

followers={'else','elseif','case','otherwise'};

finish='end';

arobase=regexp(txt,'@#\s*(?<type>\w+)\s*(?<action>.*)(?<pos>)(?<level>)','names');

for ib=1:numel(arobase)
    
    if ~isempty(arobase{ib})
        
        arobase{ib}.pos=ib;
        
    end
    
end

arobase=[arobase{:}];

if isempty(arobase)
    
    return
    
end

if any(strcmp(arobase(1).type,triggers))
    
    arobase(1).level=0;
    
else
    
    error('wrong level start')
    
end

depth=0;

for iro=2:numel(arobase)
    
    if any(strcmp(arobase(iro).type,triggers))
        
        depth=depth+1;
        
        arobase(iro).level=depth;
        
    elseif strcmp(arobase(iro).type,finish)
        
        arobase(iro).level=depth;
        
        depth=depth-1;
        
    elseif any(strcmp(arobase(iro).type,followers))
        
        arobase(iro).level=depth;
        
    end
    
end

levels=[arobase.level];

if any(levels<0)
    
    error('too many endings')
    
end

good=levels==0;

arobase=arobase(good);

types={arobase.type};

first_end=find(strcmp(types,finish),1,'first');

if isempty(first_end)
    
    error('opened flow control not closed')
    
end

arobase=arobase(1:first_end);

% check consistency
types=types(1:first_end);

theLeader=types{1};

theFollowers=types(2:end-1);

switch theLeader
    
    case 'for'
        
        if ~isempty(theFollowers)
            
            error('"For" should not have followers')
            
        end
        
    case 'if'
        
        if ~isempty(theFollowers)
            
            ellsse=find(strcmp(theFollowers,'else'));
            
            if numel(ellsse)==1
                
                if ~strcmp(types{end},'end')
                    
                    error('for "if", "else" must appear right before "end"')
                    
                end
                
                theFollowers(ellsse)=[];
                
            elseif numel(ellsse)>1
                
                error('in "if", only one "else" is allowed')
                
            end
            
            ellsseiff=strcmp(theFollowers,'elseif');
            
            theFollowers(ellsseiff)=[];
            
            if ~isempty(theFollowers)
                
                error('Followers of "if" can only be "else" and "elseif"')
                
            end
            
        end
        
    case 'switch'
        
        if isempty(theFollowers)
            
            error('"switch" must have followers')
            
        end
        
        ozerwise=find(strcmp(theFollowers,'otherwise'));
        
        if numel(ozerwise)==1
            
            if ~strcmp(theFollowers{end},'otherwise')
                
                error('for "switch", "otherwise" must appear right before "end"')
                
            end
            
            theFollowers(ozerwise)=[];
            
        elseif numel(ozerwise)>1
            
            error('in "switch", only one "otherwise" is allowed')
            
        end
        
        kases=strcmp(theFollowers,'case');
        
        theFollowers(kases)=[];
        
        if ~isempty(theFollowers)
            
            error('Followers of "switch" can only be "case" and "otherwise"')
            
        end
        
end


end

%--------------------------------------------------------------------------

%{
function xout=process_keywords(xin)

myreplace=@replacer; %#ok<NASGU>

kwords={'movave','movavg','movgeom','movprod','movsum','diff','difflog'};

kwords=parser.cell2matize(kwords);

patt=['\<',kwords,'\>\(\s*(\w+)\s*(,\s*\d+)?\s*\)'];

xout=regexprep(xin,patt,'${myreplace($1,$2,$3)}');

    function str=replacer(func,v,n)
        
        if isempty(n)
            
            n='1';
            
        end
        
        n=str2double(strrep(n,',',''));
        
        if n<=0||~isfinite(n)||(floor(n)~=ceil(n))
            
            error(['wrong specification of integer when using ',func])
            
        end
        
        switch func
            
            case {'movave','movavg'}
                
                str=many_together('+');
                
                str=['1/',int2str(n),'*',str];
                
            case 'movprod'
                
                str=many_together('*');
                
            case 'movgeom'
                
                str=many_together('*');
                
                str=[str,'^(1/',int2str(n),')'];
                
            case 'movsum'
                
                str=many_together('+');
                
            case 'diff'
                
                str=[v,'-',v,'{-',int2str(n),'}'];
                
            case 'difflog'
                
                str=['log(',v,')-log(',v,'{-',int2str(n),'})'];
                
            otherwise
                
                error('Report this error to the forum or contact junior.maih@gmail.com')
                
        end
        
        str=['(',str,')'];
        
        function str=many_together(sign_)
            
            str=['(',v];
            
            for ii=2:n
                
                str=[str,sign_,v,'{-',int2str(ii-1),'}']; %#ok<AGROW>
                
            end
            
            str=[str,')'];
            
        end
        
    end

end
%}

