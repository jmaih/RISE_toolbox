function tex_code = latex_model_file(model,model_syntax,model_line_numbers)
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

if nargin<3
    model_line_numbers=true;
    if nargin<2
        model_syntax=true;
    end
end

tex_code = '';

% par_list = prepare_list(get(model,'par_list')); % bold
% exo_list = prepare_list(get(model,'exo_list'));
% endo_list = prepare_list(get(model,'endo_list')); 
% chain_list = prepare_list(get(model,'chain_list')); % red

kwd_list=prepare_list(load_keywords()); % blue

br = sprintf('\n');

raw_file=model.raw_file;
nline = numel(raw_file);

filename='';
for ii = 1 : nline
    % Split the line if there is a line comment.
    raw_line=raw_file(ii);
    line_number=raw_line.line_numbers;
    if isempty(filename)|| ~strcmp(filename,raw_line.filename)
        filename=raw_line.filename;
        new_file=true;
    end
    tok = regexp(raw_line.code,'([^%]*)(%.*)?','tokens','once');
    if ~isempty(tok)
        x = tok{1};
        y = tok{2};
    else
        x = '';
        y = '';
    end
    if new_file
        tmp=process_line(filename);
        new_file=false;
        tex_code = [tex_code,tmp,' \\',br]; %#ok<AGROW>
    end
    x = process_line(x);
    y = '';%docomments(y);
    tex_code = [tex_code,x,y,' \\',br]; %#ok<AGROW>
end
% escape #
tex_code=strrep(tex_code,'#','\#');

    function xx = process_line(xx)
        myprescreen=@pre_screen; %#ok<NASGU>
        if model_syntax
            % keywords
            pat1=['(?<!\w)',kwd_list,'(?!\w)']; 
            replace1='@\\textcolor{blue}{\\texttt{${myprescreen($0)}}}\\verb@';
            xx = regexprep(xx,pat1,replace1);
%             % numbers
%             pat2='(?<![a-z_A-Z,\.{])(\d*\.*\d+)(?!})'; replace2='@\\textcolor{red}{\\texttt{$1}}\\verb@';
%             xx = regexprep(xx,pat2,replace2);
            % time
            pat3='(\w+)({)((+|-)*\d+)(})'; replace3='$1$2@\\textcolor{red}{\\texttt{$3}}\\$4\\verb@';
            xx = regexprep(xx,pat3,replace3);
        end
        if new_file
%             xx=['@\textbf{\texttt{',strrep(xx,'_','\_'),'}}\verb@'];
            xx=['@\textcolor{green}{\texttt{\% new file with name :: ',strrep(xx,'_','\_'),'}}\verb@'];
        end
        if model_line_numbers
            xx = [sprintf('%0.0f: ',line_number),xx];
        end
        xx = ['\verb@',xx,'@'];
    end

end
function kwd_list=load_keywords()
kwd_list={'@','#','define','dsge_var','endfor','endif','include',...
    'linear','log_vars','model','parameters','planner_discount',...
    'planner_objective','steady_state_model','var','varexo',...
    'varobs','observables','endogenous','exogenous','imposed',...
    'parameterization','parameter_restrictions','unique'...
    'exogenous_definition','else','elseif','end','for','if',...
    'switch','while'};
end

function list=prepare_list(list)
list=list(:)';
list=strcat(list,'|');
list=cell2mat(list);
list=['(',list(1:end-1),')'];
end

function c=pre_screen(c)
c=strrep(c,'_','\_');
end
