function block=capture_parameterization(dictionary,listing)
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

PARAM_DELIMITERS=[char([9:13,32]),',;'];

parameter_names={dictionary.parameters.name};
parameter_chains=[dictionary.parameters.governing_chain];
markov_chain_names={dictionary.markov_chains.name};

% initialize output
%------------------
block=struct();

% first make sure all statements are on the same line
listing=reprocess_parameter_batch(listing);
nlisting=size(listing,1);

% now do it
PP=struct();
is_estimated=false(nlisting,1);
is_dirichlet=false(nlisting,1);
error_control=cell(nlisting,2);
ncount=0;
current_names={};
push_if_validated=@parser.push_if_validated;
for ii=1:nlisting
    capture_parameterization_engine(listing(ii,:));
end

% get the baseline calibration and the priors
%--------------------------------------------
Calibration=struct();
Priors=struct();
pnames=fieldnames(PP);
for iname=1:numel(pnames)
    if is_estimated(iname)
        Priors.(pnames{iname})=PP.(pnames{iname});
    end
    if is_dirichlet(iname)
        % names in first row, values in second row
        tmp=reshape(PP.(pnames{iname})(2:end),2,[]);
        for icol=1:size(tmp,2)
            Calibration.(tmp{1,icol})=tmp{2,icol};
        end
    else
        Calibration.(pnames{iname})=PP.(pnames{iname}){1};
    end
end
% also trim the error control
%----------------------------
error_control=error_control(is_estimated,:);

% transfer everything to a single structure
%------------------------------------------
block.Calibration=Calibration;
block.Priors=Priors;
block.ncount=ncount;
block.error_control=error_control;

    function capture_parameterization_engine(cell_info)
        % the parameterization row could span several lines, in that case
        % those lines need to be put together, should an error be issued
        %-----------------------------------------------------------------
        iline_=cell_info{1};
        old_line=iline_;
        iline_=sprintf('%0.0f',old_line(1));
        for iii=2:numel(old_line)
            iline_=[iline_,' & ',sprintf('%0.0f',old_line(iii))]; %#ok<AGROW>
        end
        rawline_=cell_info{2};
        rawline_(isspace(rawline_))=[];
        % remove the semi-colon
        if strcmp(rawline_(end),';')
            rawline_=rawline_(1:end-1);
        end
        file_name_=cell_info{3};
        if isempty(rawline_)
            % I do not expect this to happen at this stage
            return
        end
        name_file_line={[],file_name_,iline_};
        error_control(ii,:)=name_file_line(2:end);
        [tokk,rawline_]=strtok(rawline_,[PARAM_DELIMITERS,'(']);
        param_location=find(strcmp(tokk,parameter_names));
        if isempty(param_location)
            is_dirichlet(ii)=strncmp('dirichlet',tokk,9);
            if ~is_dirichlet(ii)
                error([mfilename,':: ',tokk,...
                    ' is not recognized as a parameter in ',file_name_,...
                    ' at line(s) ',iline_])
            end
        end
        
        % initialize and begin storing
        %-----------------------------
        par_name=tokk;
        if strcmp(rawline_(1),'(')
            right_par=strfind(rawline_,')');
            if isempty(right_par)
                error([mfilename,':: closing parenthesis missing for parameter name specification in ',file_name_,' at line(s) ',iline_])
            end
            right_par=right_par(1);
            if isempty(rawline_(2:right_par-1))
                error([mfilename,':: content of parentheses after the parameter name should not be empty in ',file_name_,' at line(s) ',iline_])
            end
            par_name=[par_name,rawline_(1:right_par)];
            check_affiliation(par_name)
            par_name=parser.param_texname_to_param_name(par_name);
            rawline_=rawline_(right_par+1:end);
        end
        % check that the parameter is not listed twice in the same state
        %---------------------------------------------------------------
        if is_dirichlet(ii)
            old_tex_name=par_name;
            group='group';
        else
            group='';
            old_tex_name=parser.param_name_to_param_texname(par_name,markov_chain_names);
        end
        if ismember(par_name,current_names)
            error([mfilename,':: parameter ',group,' ',old_tex_name,...
                ' appears for a second time in ',file_name_,...
                ' at line(s) ',iline_])
        end
        current_names=[current_names,par_name];
        PP.(par_name)={};
        ncount=ncount+1;
        % for the rest, just count the number of commas to separate items
        %----------------------------------------------------------------
        % the semicolon was removed above
        name_file_line{1}=old_tex_name;
        if is_dirichlet(ii)
            process_dirichlet()
        else
            process_normal()
        end
        function check_affiliation(par_name)
            par_name=parser.param_name_to_param_texname(par_name,...
                markov_chain_names);
            left_par=strfind(par_name,'(');
            right_par__=strfind(par_name,')');
            style=~isempty(left_par)+2*(~isempty(right_par__));
            switch style
                case 0 % nice
                    first='const';
                    second='1';
                case {1,2}
                    error(['parenthesis mismatch for parameter name specification in ',file_name_,' at line(s) ',iline_])
                case 3 % nice
                    right_par__=right_par__(1);
                    extract=par_name(left_par+1:right_par__-1);
                    par_name=par_name(1:left_par-1);
                    % locate the comma and split
                    comma=strfind(extract,',');
                    if isempty(comma)||numel(comma)~=1
                        error([mfilename,':: exacly 2 arguments are to enter the parentheses right after the parameter name in ',file_name_,' at line(s) ',iline_])
                    end
                    first=extract(1:comma-1);
                    second=extract(comma+1:end);
            end
            chain_loc=find(strcmp(first,markov_chain_names));
            if isempty(chain_loc)
                error([' parameter ',par_name,'''s governing chain ',first,...
                    ' is not recognized as a chain name in ',file_name_,' at line(s) ',iline_])
            end
            if ~(parameter_chains(param_location)==chain_loc)
                error([' parameter ',par_name,' is governed by markov chain "',...
                    markov_chain_names{parameter_chains(param_location)},...
                    '" and not "',first,'" as in ',file_name_,...
                    ' at line(s) ',iline_])
            end
            eval_tok=eval(second);
            if ~isfinite(eval_tok)
                error([mfilename,':: parameter ',par_name,...
                    ' has nonsensical expression for its state in ',...
                    file_name_,' at line(s) ',iline_])
            end
        end
        function process_dirichlet()
            is_estimated(ii)=true;
            % make sure there is no funny parameter in the form
            % pname(chain,state). replace it with pname_chain_state
            rawline_=regexprep(rawline_,'(\w+)\((\w+),(\d+)\)','$1_$2_$3');
            odd=true;
            while ~isempty(rawline_)
                [tokk,rawline_]=strtok(rawline_,PARAM_DELIMITERS); %#ok<STTOK>
                if odd
                    dd=push_if_validated(str2double(tokk),[],'dirichlet',...
                        name_file_line);
                else
                    dd=tokk;
                    check_affiliation(tokk)
                end
                PP.(par_name)=[PP.(par_name),{dd}];
                odd=~odd;
            end
            if odd
                error(['Wrong specificiation of dirichlet in ',file_name_,...
                    ' at line(s) ',iline_])
            end
        end
        function process_normal()
            commas=strfind(rawline_,',');
            ncom=numel(commas);
            if ncom==0 % start value
                error(['could not find the start value for ',old_tex_name,...
                    ' in ',file_name_,' at line(s) ',iline_])
            else
                if ncom==1
                    start_=push_if_validated(eval(rawline_(commas(1)+1:end)),...
                        [],'start',name_file_line);
                    PP.(par_name)=[PP.(par_name),{start_}];
                else
                    if ncom<3
                        error(['an estimated parameter should have at least 3 components ',old_tex_name,...
                            ' in ',file_name_,' at line(s) ',iline_])
                    end
                    is_estimated(ii)=true;
                    start_=push_if_validated(eval(rawline_(commas(1)+1:commas(2)-1)),...
                        [],'start',name_file_line);
                    PP.(par_name)=[PP.(par_name),{start_}];
                    lowqtlOrMean=eval(rawline_(commas(2)+1:commas(3)-1));
                    lowqtlOrMean=push_if_validated(lowqtlOrMean,[],...
                        'lower quantile or mean',name_file_line);
                    PP.(par_name)=[PP.(par_name),{lowqtlOrMean}];
                    if ncom==3
                        upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:end)),...
                            [],'upper quantile or standard deviation',name_file_line);
                        PP.(par_name)=[PP.(par_name),{upperQtlOrStdev}];
                    else
                        upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:commas(4)-1)),...
                            [],'upper quantile or standard deviation',name_file_line);
                        PP.(par_name)=[PP.(par_name),{upperQtlOrStdev}];
                        if ncom==4
                            distribution_prob=rawline_(commas(4)+1:end);
                            PP.(par_name)=[PP.(par_name),{distribution_prob}];
                        else
                            distribution_prob=rawline_(commas(4)+1:commas(5)-1);
                            PP.(par_name)=[PP.(par_name),{distribution_prob}];
                            if ncom==5
                                lb=push_if_validated(eval(rawline_(commas(5)+1:end)),...
                                    [],'lower bound',name_file_line);
                                PP.(par_name)=[PP.(par_name),{lb}];
                            else
                                lb=push_if_validated(eval(rawline_(commas(5)+1:commas(6)-1)),...
                                    [],'lower bound',name_file_line);
                                PP.(par_name)=[PP.(par_name),{lb}];
                                if ncom==6
                                    ub=eval(rawline_(commas(6)+1:end));
                                    ub=push_if_validated(ub,[],'upper bound',name_file_line);
                                    PP.(par_name)=[PP.(par_name),{ub}];
                                else
                                    error(['too many commas prevent parsing in ',old_tex_name,' in ',file_name_,' at line(s) ',iline_])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function newbatch=reprocess_parameter_batch(batch)
% this is to make sure all complete statements are on the same line
newbatch=cell(0,3);
lines=[];
dough='';
while ~isempty(batch)
    bb=batch(1,:);
    lines=[lines,bb{1}]; %#ok<AGROW>
    dough=[dough,bb{2}]; %#ok<AGROW>
    dough=strtrim(dough);
    if strcmp(dough(end),';')
        newbatch=[newbatch;
            {lines,dough,bb{3}}]; %#ok<AGROW>
        lines=[];
        dough='';
    end
    batch=batch(2:end,:);
end
end