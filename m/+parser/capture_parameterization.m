function block=capture_parameterization(dictionary,listing)
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
error_control=cell(nlisting,2);
ncount=0;
for ii=1:nlisting
    capture_parameterization_engine(listing(ii,:));
end

% get the baseline calibration and the priors
%--------------------------------------------
Calibration=struct();
Priors=struct();
pnames=fieldnames(PP);
for iname=1:numel(pnames)
    Calibration.(pnames{iname})=PP.(pnames{iname}){1};
    if is_estimated(iname)
        Priors.(pnames{iname})=PP.(pnames{iname});
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
        error_control(ii,:)={file_name_,iline_};
        [tokk,rawline_]=strtok(rawline_,[PARAM_DELIMITERS,'(']);
        param_location=find(strcmp(tokk,parameter_names));
        if isempty(param_location)
            error([mfilename,':: ',tokk,' is not recognized as a parameter in ',file_name_,' at line(s) ',iline_])
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
            extract=rawline_(2:right_par-1);
            if ~isempty(extract)
                % locate the comma and split
                comma=strfind(extract,',');
                if isempty(comma)||numel(comma)~=1
                    error([mfilename,':: exacly 2 arguments are to enter the parentheses right after the parameter name in ',file_name_,' at line(s) ',iline_])
                else
                    first=extract(1:comma-1);
                    second=extract(comma+1:end);
                    chain_loc=find(strcmp(first,markov_chain_names));
                    if isempty(chain_loc)
                        error([' parameter ',par_name,'''s governing chain ',first,...
                            ' is not recognized as a chain name in ',file_name_,' at line(s) ',iline_])
                    end
                    if ~(parameter_chains(param_location)==chain_loc)
                        error([' parameter ',par_name,' is governed by markov chain "',markov_chain_names(parameter_chains(param_location)),...
                            '" and not "',first,'" as in ',file_name_,' at line(s) ',iline_])
                    end
                        eval_tok=eval(second);
                    if ~isfinite(eval_tok)
                        error([mfilename,':: parameter ',par_name,' has nonsensical expression for its state in ',file_name_,' at line(s) ',iline_])
                    end
                end
            else
                error([mfilename,':: content of parentheses after the parameter name should not be empty in ',file_name_,' at line(s) ',iline_])
            end
            par_name=[par_name,rawline_(1:right_par)];
            par_name=parser.param_name_to_valid_param_name(par_name);
            rawline_=rawline_(right_par+1:end);
        end
        % check that the parameter is not listed twice in the same state
        %---------------------------------------------------------------
        old_tex_name=parser.valid_param_name_to_tex_name(par_name,markov_chain_names);
        if ncount
            current_names=fieldnames(PP);
            if ismember(par_name,current_names)
                error([mfilename,':: parameter ',old_tex_name,' appears for a second time in ',file_name_,' at line(s) ',iline_])
            end
        end
        PP.(par_name)={};
        ncount=ncount+1;
        % for the rest, just count the number of commas to separate items
        %----------------------------------------------------------------
        % the semicolon was removed above
        commas=strfind(rawline_,',');
        ncom=numel(commas);
        if ncom==0 % start value
            error(['could not find the start value for ',old_tex_name,...
                ' in ',file_name_,' at line(s) ',iline_])
        else
            if ncom==1
                start_=push_if_validated(eval(rawline_(commas(1)+1:end)),'start');
                PP.(par_name)=[PP.(par_name),{start_}];
            else
                is_estimated(ii)=true;
                start_=push_if_validated(eval(rawline_(commas(1)+1:commas(2)-1)),'start');
                PP.(par_name)=[PP.(par_name),{start_}];
                lowqtlOrMean=eval(rawline_(commas(2)+1:commas(3)-1));
                lowqtlOrMean=push_if_validated(lowqtlOrMean,'lower quantile or mean');
                PP.(par_name)=[PP.(par_name),{lowqtlOrMean}];
                if ncom==3
                    upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:end)),'upper quantile or standard deviation');
                    PP.(par_name)=[PP.(par_name),{upperQtlOrStdev}];
                else
                    upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:commas(4)-1)),'upper quantile or standard deviation');
                    PP.(par_name)=[PP.(par_name),{upperQtlOrStdev}];
                    if ncom==4
                        distribution_prob=rawline_(commas(4)+1:end);
                        PP.(par_name)=[PP.(par_name),{distribution_prob}];
                    else
                        distribution_prob=rawline_(commas(4)+1:commas(5)-1);
                        PP.(par_name)=[PP.(par_name),{distribution_prob}];
                        if ncom==5
                            lb=push_if_validated(eval(rawline_(commas(5)+1:end)),'lower bound');
                            PP.(par_name)=[PP.(par_name),{lb}];
                        else
                            lb=push_if_validated(eval(rawline_(commas(5)+1:commas(6)-1)),'lower bound');
                            PP.(par_name)=[PP.(par_name),{lb}];
                            if ncom==6
                                ub=eval(rawline_(commas(6)+1:end));
                                ub=push_if_validated(ub,'upper bound');
                                PP.(par_name)=[PP.(par_name),{ub}];
                            else
                                error(['too many commas prevent parsing in ',old_tex_name,' in ',file_name_,' at line(s) ',iline_])
                            end
                        end
                    end
                end
            end
        end
        function dd=push_if_validated(val,type)
            if isnan(val)
                error(['wrong specification of ',type,' value for ',old_tex_name,' in ',file_name_,' at line(s) ',iline_])
            else
                dd=val;
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

% % function block=capture_parameterization(dictionary,listing)
% % PARAM_DELIMITERS=[char([9:13,32]),',;'];
% % 
% % parameter_names={dictionary.parameters.name};
% % parameter_chains=[dictionary.parameters.governing_chain];
% % markov_chain_names={dictionary.markov_chains.name};
% % % first make sure all statements are on the same line
% % listing=reprocess_parameter_batch(listing);
% % % now do it
% % block=[];
% % for ii=1:size(listing,1)
% %     block=capture_parameterization_engine(block,listing(ii,:));
% % end
% % 
% % % check that every parameter is controled by one chain only and that it is
% % % assigned a value in every state of the commanding chain. NB: above, I
% % % have already made sure that every parameter in the parameterization block
% % % is declared
% % % we set the default governing chain to const rather than nan
% % for ii=1:numel(parameter_names)
% %     par_i=parameter_names{ii};
% %     loc=find(strcmp(par_i,{block.name}));
% %     if ~isempty(loc)
% %         % check that there no duplicate statements
% %         for j1=1:numel(loc)
% %             for j2=j1+1:numel(loc)
% %                 if isequal(block(loc(j1)).state,block(loc(j2)).state)
% %                     error([mfilename,':: parameterization of ',par_i,' occurs at least twice in the same state'])
% %                 end
% %             end
% %         end
% %     end
% % end
% % 
% %     function block=capture_parameterization_engine(block,cell_info)
% %         if nargin==0
% %             block=struct('name','','chain','const','state',1,'start',nan,...
% %                 'lower_quantile',nan,'upper_quantile',nan,'prior_mean',nan,...
% %                 'prior_stdev',nan,'prior_distrib','uniform','prior_prob',1,...
% %                 'lower_bound',nan,'upper_bound',nan);
% %             return
% %         end
% %         if isempty(block)
% %             block=struct('name',{},'chain',{},'state',{},'start',{},...
% %                 'lower_quantile',{},'upper_quantile',{},'prior_mean',{},...
% %                 'prior_stdev',{},'prior_distrib',{},'prior_prob',{},...
% %                 'lower_bound',{},'upper_bound',{});
% %         end
% %         iline_=cell_info{1};
% %         old_line=iline_;
% %         iline_=sprintf('%0.0f',old_line(1));
% %         for iii=2:numel(old_line)
% %             iline_=[iline_,' & ',sprintf('%0.0f',old_line(iii))];
% %         end
% %         rawline_=cell_info{2};
% %         rawline_(isspace(rawline_))=[];
% %         % remove the semi-colon
% %         if strcmp(rawline_(end),';')
% %             rawline_=rawline_(1:end-1);
% %         end
% %         file_name_=cell_info{3};
% %         if isempty(rawline_)
% %             return
% %         end
% %         [tokk,rawline_]=strtok(rawline_,[PARAM_DELIMITERS,'(']);
% %         param_location=find(strcmp(tokk,parameter_names));
% %         if isempty(param_location)
% %             error([mfilename,':: ',tokk,' is not recognized as a parameter in ',file_name_,' at line(s) ',iline_])
% %         end
% %         
% %         % initialize and begin storing
% %         %-----------------------------
% %         %         nblks=numel(block)+1;
% %         newblk=capture_parameterization_engine();
% %         newblk.name=tokk;
% %         if strcmp(rawline_(1),'(')
% %             right_par=strfind(rawline_,')');
% %             right_par=right_par(1);
% %             extract=rawline_(2:right_par-1);
% %             rawline_=rawline_(right_par+1:end);
% %             if ~isempty(extract)
% %                 % locate the comma and split
% %                 comma=strfind(extract,',');
% %                 if isempty(comma)||numel(comma)~=1
% %                     error([mfilename,':: exacly 2 arguments are to enter the parentheses right after the parameter name in ',file_name_,' at line(s) ',iline_])
% %                 else
% %                     first=extract(1:comma-1);
% %                     second=extract(comma+1:end);
% %                     chain_loc=find(strcmp(first,markov_chain_names));
% %                     if isempty(chain_loc)
% %                         error([' parameter ',newblk.name,'''s governing chain ',first,...
% %                             ' is not recognized as a chain name in ',file_name_,' at line(s) ',iline_])
% %                     end
% %                     newblk.chain=first;
% %                     if ~(parameter_chains(param_location)==chain_loc)
% %                         error([' parameter ',newblk.name,' is governed by markov chain "',markov_chain_names(parameter_chains(param_location)),...
% %                             '" and not "',first,'" as in ',file_name_,' at line(s) ',iline_])
% %                     end
% %                     try
% %                         eval_tok=eval(second);
% %                     catch %#ok<*CTCH>
% %                         error([mfilename,':: parameter ',newblk.name,' has nonsensical expression for its state in ',file_name_,' at line(s) ',iline_])
% %                     end
% %                     newblk.state=eval_tok;
% %                 end
% %             else
% %                 error([mfilename,':: content of parentheses after the parameter name should not be empty in ',file_name_,' at line(s) ',iline_])
% %             end
% %         end
% %         % for the rest, just count the number of commas to separate items
% %         %----------------------------------------------------------------
% %         % the semicolon was removed above
% %         commas=strfind(rawline_,',');
% %         ncom=numel(commas);
% %         if ncom==0 % start value
% %             error(['could not find the start value for ',newblk.name,' in ',file_name_,' at line(s) ',iline_])
% %         else
% %             if ncom==1
% %                 newblk.start=push_if_validated(eval(rawline_(commas(1)+1:end)),'start');
% %             else
% %                 newblk.start=push_if_validated(eval(rawline_(commas(1)+1:commas(2)-1)),'start');
% %                 if ncom==2
% %                     error(['when the lower quantile(mean) is specified, the upper quantile(stdev) must be specified. for ',newblk.name,' in ',file_name_,' at line(s) ',iline_])
% %                 end
% %                 lowqtlOrMean=eval(rawline_(commas(2)+1:commas(3)-1));
% %                 lowqtlOrMean=push_if_validated(lowqtlOrMean,'lower quantile or mean');
% %                 if ncom==3
% %                     upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:end)),'upper quantile or standard deviation');
% %                     % push the above into the quantiles for a uniform
% %                     % distribution. We keep the probability to 1
% %                     newblk.lower_quantile=lowqtlOrMean;
% %                     newblk.upper_quantile=upperQtlOrStdev;
% %                 else
% %                     upperQtlOrStdev=push_if_validated(eval(rawline_(commas(3)+1:commas(4)-1)),'upper quantile or standard deviation');
% %                     if ncom==4
% %                         distribution_prob=rawline_(commas(4)+1:end);
% %                     else
% %                         distribution_prob=rawline_(commas(4)+1:commas(5)-1);
% %                         if ncom==5
% %                             lb=eval(rawline_(commas(5)+1:end));
% %                         else
% %                             lb=eval(rawline_(commas(5)+1:commas(6)-1));
% %                             if ncom==6
% %                                 ub=eval(rawline_(commas(6)+1:end));
% %                                 newblk.upper_bound=push_if_validated(ub,'upper bound');
% %                             else
% %                                 error(['too many commas prevent parsing in ',newblk.name,' in ',file_name_,' at line(s) ',iline_])
% %                             end
% %                         end
% %                         newblk.lower_bound=push_if_validated(lb,'lower bound');
% %                     end
% %                     left_par=strfind(distribution_prob,'(');
% %                     if isempty(left_par)
% %                         left_par=length(distribution_prob)+1;
% %                         newblk.prior_mean=lowqtlOrMean;
% %                         newblk.prior_stdev=upperQtlOrStdev;
% %                         % we have to change default for the probability
% %                         newblk.prior_prob=nan;
% %                     else
% %                         right_par=strfind(distribution_prob,')');
% %                         the_prob=eval(distribution_prob(left_par+1:right_par-1));
% %                         the_prob=push_if_validated(the_prob,'probability');
% %                         if the_prob<0||the_prob>1
% %                             error(['probability value should be between 0 and 1 for ',newblk.name,' in ',file_name_,' at line(s) ',iline_])
% %                         end
% %                         newblk.prior_prob=the_prob;
% %                         newblk.lower_quantile=lowqtlOrMean;
% %                         newblk.upper_quantile=upperQtlOrStdev;
% %                     end
% %                     newblk.prior_distrib=distribution_prob(1:left_par-1);
% %                 end
% %             end
% %         end
% %         block(end+1)=newblk;
% %         function dd=push_if_validated(val,type)
% %             if isnan(val)
% %                 error(['wrong specification of ',type,' value for ',newblk.name,' in ',file_name_,' at line(s) ',iline_])
% %             else
% %                 dd=val;
% %             end
% %         end
% %     end
% % 
% % end
% % 
% % function newbatch=reprocess_parameter_batch(batch)
% % % this is to make sure all complete statements are on the same line
% % newbatch=cell(0,3);
% % lines=[];
% % dough='';
% % while ~isempty(batch)
% %     bb=batch(1,:);
% %     lines=[lines,bb{1}];
% %     dough=[dough,bb{2}];
% %     dough=strtrim(dough);
% %     if strcmp(dough(end),';')
% %         newbatch=[newbatch;
% %             {lines,dough,bb{3}}];
% %         lines=[];
% %         dough='';
% %     end
% %     batch=batch(2:end,:);
% % end
% % end
