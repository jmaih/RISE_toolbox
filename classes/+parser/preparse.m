function [output,has_macro]=preparse(FileName,definitions)
% int2str(x)=sprintf('%.0f',x)
if nargin<2
    definitions=struct();
end

valid_extensions={'rs','dyn','mod'};
FileName(isspace(FileName))=[];
loc=strfind(FileName,'.');
if isempty(loc)
    for iext=1:numel(valid_extensions)
        found= exist([FileName,'.',valid_extensions{iext}],'file');
        if found
            break
        end
    end
else
    ext=FileName(loc+1:end);
    if ~ismember(ext,valid_extensions)
        error(['file extension ',ext,' is not recognized by RISE'])
    end
    found= exist(FileName,'file');
end
if ~found
    error([mfilename,':: ',FileName,'.rs or ',FileName,'.dyn or ',FileName,'.mod not found'])
end

% first sweep: read the file, and remove comments and stamp the lines
RawFile=parser.read_file(FileName);

% second sweep: process the if else end for and include comments
% include should also go through preparsing
[output,has_macro]=preparse_expand(RawFile,definitions);

end

function [output,has_macro]=preparse_expand(RawFile,definitions,has_macro,output)
if nargin<4
    output=[];
    if nargin<3
        has_macro=false;
    end
end
if isempty(RawFile)
    return
end
if isempty(output),output=cell(0,3);end
if isempty(has_macro),has_macro=false;end

SPACE_DELIMITERS=char([9:13,32]);

depth=0;
iter=0;
end_of_file=size(RawFile,1);
while iter<end_of_file
    iter=iter+1;
    rline=next_line(iter);
    switch rline.first_tok
        case '@#include'
            has_macro=true;
            if depth==0
                output=parser.append_file(output,rline,definitions);
            end
        case {'@#if','@#for'}
            has_macro=true;
            depth=depth+1;
            control=rline.first_tok(3:end);
            if depth==1
                if strcmp(control,'if')
                    steps=check_and_store_condition(rline,definitions,iter);
                    % steps={start,finish,flag};
                else
                    steps=get_index_start_and_finish(rline,definitions,iter);
                    % steps={start,finish,{index,start_loop,end_loop}};
                end
                jter=iter;
                old_line=rline;
                while jter<end_of_file
                    jter=jter+1;
                    rline=next_line(jter);
                    update_steps_and_depth(rline);
                    if depth==0
                        break
                    end
                end
                if depth~=0
                    error(['if statement started in ',old_line.FileName,' at line ',old_line.row_number_string,' is not closed'])
                end
                step_rows=select_appropriate_rows(steps);
                if strcmp(control,'if')
                    [output,has_macro]=preparse_expand(RawFile(step_rows,:),definitions,has_macro,output);
                else
                    [output,has_macro2]=for_loop_batch(output,RawFile(step_rows,:),definitions,steps{3});
                    has_macro=has_macro||has_macro2;
                end
                iter=jter;
            end
        case {'@#elseif','@#else','@#end'}
            error(['token ',rline.first_tok,' in ',rline.FileName,' at line ',rline.row_number_string,' is not initialized'])
        otherwise
            tank=parser.remove_definitions(RawFile(iter,:),definitions);
            output=[output;tank];
    end
end

    function update_steps_and_depth(rline)
        switch rline.first_tok
            case {'@#if','@#for'}
                depth=depth+1;
            case '@#elseif'
                if depth==1
                    newsteps=check_and_store_condition(rline,definitions,jter);
                    steps{end,2}=jter-1;
                    steps=[steps;newsteps];
                    good=find(cell2mat(steps(:,end)));
                    if numel(good)>1
                        error(['multiple valid conditional statements in "if" control in ',rline.FileName,' at line ',rline.row_number_string])
                    end
                end
            case '@#else'
                if depth==1
                    steps{end,2}=jter-1;
                    flag=~any(cell2mat(steps(:,end)));
                    steps=[steps;{jter+1,nan,flag}];
                end
            case {'@#end','@#endif','@#endfor'} % accommodating dynare
                if depth==1
                    steps{end,2}=jter-1;
                end
                depth=depth-1;
            otherwise
        end
    end

%--------------------------------------------------------------------------
    function rline=next_line(iter)
        rline.rawline = RawFile{iter,1};
        rline.FileName= RawFile{iter,2};
        rline.row_number_string= int2str(RawFile{iter,3});
        
        [tokk,rest_]=strtok(rline.rawline ,SPACE_DELIMITERS);
        if strncmp(tokk,'@#',2)
            if length(tokk)==2
                oldtokk=tokk;
                [tokk,rest_]=strtok(rest_,SPACE_DELIMITERS);
                tokk=[oldtokk,tokk];
            end
        end
        rline.first_tok=tokk;
        rline.rest=rest_;
    end
end

%--------------------------------------------------------------------------
function step_rows=select_appropriate_rows(steps)
step_rows=1:0;
if iscell(steps{end,end})
    % we have a for loop
    step_rows=steps{1}:steps{2};
else
    % we have an if control
    good=find(cell2mat(steps(:,end)));
    if ~isempty(good)
        steps=steps(good,:);
        step_rows=steps{1}:steps{2};
    end
end
end
%--------------------------------------------------------------------------
function [output,has_macro]=for_loop_batch(output,batch,definitions,index_start_finish)

index=index_start_finish{1};
start_finish=index_start_finish{2};
co=index(3:end-1);
if ischar(start_finish)
    start_finish=definitions.(start_finish);
elseif isa(start_finish,'double')
    start=start_finish(1);
    finish=start_finish(2);
    start_finish=num2cell(start:finish);
else
    error(['unknown case ',start_finish])
end
for ii=1:numel(start_finish)
    if isa(start_finish{ii},'double')
        start_finish{ii}=num2str(start_finish{ii});
    end
    definitions.(co)=start_finish{ii};
    %------------------------------------
    tank=parser.remove_definitions(batch,definitions);
    %------------------------------------
    % Now expand the content of the loop before storing...
    [tank,has_macro]=preparse_expand(tank,definitions);
    output=[output;tank];
end

end
%--------------------------------------------------------------------------
function steps=get_index_start_and_finish(rline,definitions,iter)

SPACE_DELIMITERS=char([9:13,32]);
% load definitions
fields=fieldnames(definitions);
for ifield=1:numel(fields)
    ff=fields{ifield};
    tmp=definitions.(ff); %#ok<NASGU>
    eval([ff,'=tmp;'])
end

[index,rline.rest]=strtok(rline.rest,SPACE_DELIMITERS);
index=['@{',index,'}'];
in=strfind(rline.rest,'in');
rline.rest=rline.rest(in+2:end);
colon=strfind(rline.rest,':');
if isempty(colon)
    % then we probably have a list... why not input it as {'A','adscs'}?
    start_finish=strtok(rline.rest,SPACE_DELIMITERS);
else
    start=eval(rline.rest(1:colon-1));
    finish=eval(rline.rest(colon+1:end));
    start_finish=[start,finish];
    if any(isnan(start_finish))
        error(['problem parsing ',rline.rest,' in ',rline.FileName,' at line ',rline.row_number_string])
    end
end
steps={iter+1,nan,{index,start_finish}};

end
%--------------------------------------------------------------------------
function steps=check_and_store_condition(rline,definitions,iter)
% steps={start,finish,flag};
fields=fieldnames(definitions);
for ifield=1:numel(fields)
    ff=fields{ifield};
    tmp=definitions.(ff); %#ok<NASGU>
    eval([ff,'=tmp;'])
end
flag=logical(eval(rline.rest));
if isnan(flag)
    error(['expression ',rline.rest,...
        ' in ',rline.rline.FileName,' at line ',rline.rline.row_number_string,...
        ' is not a valid logical statement'])
end
start=iter+1;
finish=nan;
steps={start,finish,flag};
% (rline.rest,definitions,iter,,)
end

