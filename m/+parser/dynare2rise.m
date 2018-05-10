function paramFileName=dynare2rise(dynFileName,riseFileName,stderr_name,do_paramFile)

if nargin<4
    
    do_paramFile=[];
    
    if nargin<3
        
        stderr_name=[];
        
        if nargin<2
            
            riseFileName=[];
            
        end
        
    end
    
end

paramFileName='';

if isempty(do_paramFile)
    
    do_paramFile=true;
    
end


if isempty(riseFileName)
    
    riseFileName=[regexprep(dynFileName,'(\w+)\.\w+','$1'),'.rs'];
    
end

if isempty(stderr_name)
    
    stderr_name='std';
    
end

dynFileName = strtrim(dynFileName);

riseFileName = strtrim(riseFileName);

% read input file
%-----------------
raw_code = read_file(dynFileName);

% Convert char(10) to white space
%--------------------------------
eol=char(10);% sprintf('\n')
% tab=char(9);% sprintf('\t')

% insert all subfiles
%--------------------
insert_all_subfiles();

% remove block comments
%-----------------------
rise_code=stepwise_removal_of_block_comments(raw_code);
% rise_code = regexprep(raw_code,'/\*(.*)\*/',''); % to revisit
lraw=length(raw_code);
lrise=length(rise_code);
fprintf(1,'Removing block comments: Old(%0.0f), New(%0.0f)\n',lraw,lrise);

header = sprintf('Conversion of Dynare file [%s] into RISE file [%s]\n',...
    dynFileName,riseFileName);

% Remove line comments.
%----------------------
lrise_old=lrise;
rise_code = regexprep(rise_code,'(//|%)(.*?\n)','\n');
lrise=length(rise_code);
fprintf(1,'Removing comment lines: Old(%0.0f), New(%0.0f)\n',lrise_old,lrise);

% replace E-015 with 10^15
%-------------------------
rise_code = regexprep(rise_code,'(\d+)(E-0)','$1/10^');

% replace y ${y}$ (long_name='output') with y "{y}(output)"
%----------------------------------------------------------
rise_code=replace_descriptions(rise_code);

% place a ¤ at the end of each @#...
%-------------------------------------------
atPound='@\s*#\s*[^\n]+';
express=['(',atPound,')\n'];
rise_code = regexprep(rise_code,express,'¤$1¤');

% replace all $ with "
%---------------------
rise_code=strrep(rise_code,'$','"');

% replace all endif or endfor with end
%-------------------------------------
rise_code=regexprep(rise_code,'(#\s*)\<end(if|for)\>','$1end');

% add space to all # signs
%----------------------------
rise_code=strrep(rise_code,'#','# ');

% remove multiple white characters
%----------------------------------
rise_code = regexprep(rise_code,'\s+',' ');

% predetermined variables
%-------------------------
[~,pred_vars]=extract_declaration_block('predetermined_variables');

% extract block of endogenous
%-----------------------------
is_endo_decl=true;
[endo_block]=extract_declaration_block('var','endogenous',is_endo_decl);

% extract block of exogenous
%-----------------------------
[exo_block,exo_list]=extract_declaration_block('varexo','exogenous');

% extract block of parameters
%-----------------------------
[param_block,par_list]=extract_declaration_block('parameters');

% extract block of observables
%-----------------------------
[obs_block,obs_list]=extract_declaration_block('varobs','observables');

% extract planner objective
%---------------------------
[planner_block]=extract_planner_objective();

% add new parameters at the end of the parameters block. Those parameters
% will not have definitions. But what if the block is empty?
%-----------------------------------------------------------------------
newparams=cell2mat(strcat([stderr_name,'_'],exo_list,'#'));

newparams=strrep(newparams,'#',' ');

param_block=[param_block,' ',newparams];

% model equations
%-----------------
rise_code = regexprep(rise_code,'model\(linear|use_dll\)','model');

% replace weird numbers that RISE does not parse yet
%---------------------------------------------------
rise_code = regexprep(rise_code,'(\d+)?(\.\d+)?e(\+|\-)(\d+)','$1$2*10^($3$4)');

rise_code = regexprep(rise_code,'(\+\s*\-|\-\s*\+)','-');

do_replace_time=false;
model_eqtns = extract_other_blocks('model;',do_replace_time);
% comment out equation tags and push the equation to the next line
%-----------------------------------------------------------------
model_eqtns=regexprep(model_eqtns,'(\[name=[^\]]*\])\s*','% $1\n');
% replace STEADY_STATE
model_eqtns=regexprep(model_eqtns,'STEADY_STATE','steady_state');
% replace ln with log
model_eqtns=regexprep(model_eqtns,'\<ln\(','log(');

% take care of predetermined variables
do_predetermined();

% steady state model equations
%------------------------------
ssmodel_eqtns = extract_other_blocks('(steady_state_model;|initval;)');

if ~isempty(ssmodel_eqtns)
    
    ssmodel_eqtns=strrep(ssmodel_eqtns,'initval;','');
    
    ssmodel_eqtns=strrep(ssmodel_eqtns,'steady_state_model;','');
    
end

% push standard deviations into equations directly
%-------------------------------------------------
add_stdev_to_model()

% Looking for covariances: shocks section
%----------------------------------------
shocks_block = extract_other_blocks('shocks;');

% Looking for covariances: shocks section
%----------------------------------------
estim_block = extract_other_blocks('estimated_params;');

est=extract_estimation(estim_block,stderr_name,obs_list);

match = regexpi(shocks_block,'corr\s*\w+\s*,\s*\w+\s*=.*?;','match');

for ii = 1 : length(match)
    
    header = [header,sprintf('Cross-covariance ignored:\n  %s\n',...
        match{ii})]; %#ok<AGROW>
    
end

% Create and save RISE code.
%---------------------------
timestamp = datestr(now);

header = [header,sprintf('\n Done %s.',timestamp)];

recreate_code();

write_parameter_file()

    function raw_code=stepwise_removal_of_block_comments(raw_code)
        
        starts=strfind(raw_code,'/*'); %starts=regexp(raw_code,'^/\*','start'); % starts=regexp(raw_code,'(?<!/)/\*','start'); 
        
        finishes=strfind(raw_code,'*/');%finishes=regexp(raw_code,'^\*/','start'); % 
        
        if numel(starts)~=numel(finishes)
            
            error('mismatches of block comments')
            
        end
        
        while ~ isempty(starts)
            
            blk=starts(end):finishes(end)+1;
            
            raw_code(blk)=[];
            
            starts(end)=[];
            
            finishes(end)=[];
            
        end
        
    end

    function write_parameter_file()
        
        if ~do_paramFile
            
            return
            
        end
        
        new_par_list=regexp(newparams,'\w+','match');
        
        par_list=[par_list,new_par_list];
        
        % parameter values
        %------------------
        % express=['\<',parser.cell2matize(par_list),'\>\s*=[^;]+;'];
        % str=regexp(rise_code,express,'match');
        
        express=['(?<name>\<',parser.cell2matize(par_list),'\>)\s*=(?<value>[^;]+);'];
        
        % first delete all places where equalities can occur i.e. model and
        % steady_state_model blocks
        tmp=regexprep(rise_code,'(model|steady_state_model)\s*(.*?)end;','');
        
        str=regexp(tmp,express,'names');
        
        pnames={str.name};
        
        pvals={str.value};
        
        % add shock standard deviations
        shock_standard_deviations();
        
        % taking care of recursive computations
        xpress=['\<',parser.cell2matize(pnames),'\>'];
        
        repl='p.$1';
        
        pvals=regexprep(pvals,xpress,repl);
        
        pvals=cellfun(@(x)x(~isspace(x)),pvals,'uniformOutput',false);
        
        code=strcat('p.',pnames,'=',pvals,';',sprintf('\n'));
        
        paramFileName=regexprep(riseFileName,'(\w+)\.\w+','$1_params');
        
        manualAdjust=[regexprep(['% ',header],'(\n)','$1%'),sprintf('\n'),...
            '% ',sprintf('\n'),...
            '% Remarks: ',sprintf('\n'),...
            '% - The parameters and shocks variances or ',sprintf('\n'),...
            '%   standard deviations found in the dynare file are assigned. ',...
            sprintf('\n'),...
            '% - Shock standard deviations not assigned in the ',sprintf('\n'),...
            '%   dynare file get a value of 0 following dynare''s convention',...
            sprintf('\n'),...
            '% - RISE will set all other parameters without a value to nan'
            ];
        
        codest=['priors=struct();',eol,eol];
        
        for ip=1:numel(est)
            
            par_i=est(ip);
            
            if strcmp(par_i.distr,'uniform')
                
                newline=sprintf('priors.%s={%s, %s, %s',...
                    par_i.name,par_i.start,par_i.lb,par_i.ub);
                
            else
                
                newline=sprintf('priors.%s={%s, %s, %s, ''%s''',...
                    par_i.name,par_i.start,par_i.mean,par_i.sd,par_i.distr);
                
                if ~isempty(par_i.lb)
                    
                    newline=sprintf('%s, %s',newline,par_i.lb);
                    
                    if ~isempty(par_i.ub)
                        
                        newline=sprintf('%s, %s',newline,par_i.ub);
                        
                    end
                    
                end
                
            end
            
            codest=[codest,newline,'};',eol,eol]; %#ok<AGROW>
            
        end
        
        code=[sprintf('function [p,priors]=%s()\n%s\n\np=struct();\n',...
            paramFileName,manualAdjust),...
            code,eol,eol,codest];
        
        %         code = [
        %             regexprep(['%',header],'(\n)','$1%'),...
        %             code
        %             ];
        
        parser.write2file(code,[paramFileName,'.m'])
        
        function shock_standard_deviations()
            % match variances
            pat=['\<var\>\s+(?<shock>',parser.cell2matize(exo_list),...
                ')\s*=(?<variance>[^;]+);'];
            
            match_shock_variances = regexp(shocks_block,pat,'names');
            
            shock_names={};
            
            shock_stdev={};
            
            if ~isempty(match_shock_variances)
                
                for ishock=1:numel(match_shock_variances)
                    
                    match_shock_variances(ishock).stdev=...
                        ['sqrt(',match_shock_variances(ishock).variance,')'];
                    
                end
                
                match_shock_stdev0=rmfield(match_shock_variances,'variance');
                
                shock_names=[shock_names,{match_shock_stdev0.shock}];
                
                shock_stdev=[shock_stdev,{match_shock_stdev0.stdev}];
                
            end
            
            % match standard deviations
            pat=['\<var\>\s+(?<shock>',parser.cell2matize(exo_list),...
                ')\s*;\s*stderr\s+(?<stdev>[^;]+);'];
            
            match_shock_stdev = regexp(shocks_block,pat,'names');
            
            if ~isempty(match_shock_stdev)
                
                shock_names=[shock_names,{match_shock_stdev.shock}];
                
                shock_stdev=[shock_stdev,{match_shock_stdev.stdev}];
                
            end
            
            if ~isempty(shock_names)
                
                shock_params=strcat([stderr_name,'_'],shock_names);
                
                pnames=[pnames,shock_params];
                
                pvals=[pvals,shock_stdev];
                
            end
            
            % set to zero the shocks that have not been assigned
            not_assigned_shocks=new_par_list(~ismember(new_par_list,pnames));
            pnames=[pnames,not_assigned_shocks];
            pvals=[pvals,repmat({'0'},1,numel(not_assigned_shocks))];
            
        end
        
    end

    function [planner_block,isCommitment,isDiscretion,discount]=extract_planner_objective()
        
        planner_block='';
        
        isDiscretion=false;
        
        isCommitment=false;
        
        discount='';
        
        loc=strfind(rise_code,'planner_objective');
        
        if isempty(loc)
            
            return
            
        end
        
        start=loc(1)+length('planner_objective');
        
        finish=find(rise_code(start:end)==';',1,'first');
        
        planner_block=rise_code(start:start+finish-1);
        
        loc=strfind(rise_code,'planner_discount');
        
        if isempty(loc)
            
            discount='0.99';
            
        else
            
            start=loc(1)+length('planner_discount');
            
            start=start+find(rise_code(start:end)=='=',1,'first');
            
            discount=strtok(rise_code(start:end),parser.delimiters);
            
        end
        
        isCommitment=~isempty(strfind(rise_code,'ramsey_policy'));
        
        isDiscretion=false;
        
        if ~isCommitment
            
            isDiscretion=~isempty(strfind(rise_code,'discretionary_policy'));
            
        end
        
        add_on=['discount=',discount,'} '];
        
        if isCommitment
            
            add_on=['{commitment=1,',add_on];
            
        elseif isDiscretion
            
            add_on=['{commitment=0,',add_on];
            
        else
            
            add_on=['{',add_on];
            
        end
        
        planner_block=['planner_objective ',add_on,planner_block];
        
    end

    function [blk,list]=extract_declaration_block(trigger,repl_trigger,is_endo_decl)
        
        if nargin<3
            % var can appear both in the declaration of endogenous and in
            % the shocks block            
            is_endo_decl=false;
            
        end
        
        express_=['(\<',trigger,'\>\s+[^;]+;)'];
        
        pull_from=rise_code;
        
        if is_endo_decl
            
            pull_from=regexprep(pull_from,'\<shocks\>;\s*(.*?)end;','');
            
        end
        
        blk=regexp(pull_from,express_,'tokens');
        
        if isempty(blk)
            
            blk='';
            
            list={};
            
            return
            
        end
        
        blk=[blk{:}];
        
        if nargin>1
            
            blk=regexprep(blk,['\<',trigger,'\>'],repl_trigger);
            
        else
            
            repl_trigger=trigger;
            
        end
        
        for iblk=1:numel(blk)
            
            semcol=find(blk{iblk}==';',1,'last');
            
            blk{iblk}(semcol)=[];
            
            if iblk>1
                blk{iblk}=strrep(blk{iblk},repl_trigger,'');
            end
            
        end
        
        blk=cell2mat(blk);% blk=blk{1};
        
        if nargout>1
            
            list=extract_list_of_names();
            
        end
        
        blk=remove_declaration_anchors(blk);
        
        function x=remove_declaration_anchors(x)
            
            x=strrep(x,'¤',sprintf('\n'));
            
        end
        
        function list=extract_list_of_names()
            
            tmp=blk;
            
            description_removal()
            
            % remove trigger
            %----------------
            loc=strfind(tmp,repl_trigger);
            
            tmp(loc(1):loc(1)+length(repl_trigger)-1)=[];
            
            % remove ¤@#...¤
            %---------------
            
            tmp=regexprep(tmp,'¤\s*@\s*#\s*[^¤\n]+¤','');
            
            list=regexp(tmp,'\<\w+\>','match'); % list=regexp(tmp,'\w+[^\s]*','match');
            
            function description_removal()
                
                dblq=find(tmp=='"');
                
                kill_spots(dblq)
                
            end
            
            function kill_spots(dblq)
                
                if ~isempty(dblq)
                    
                    dblq=reshape(dblq,2,[]);
                    
                    for iii=size(dblq,2):-1:1
                        
                        tmp(dblq(1,iii):dblq(2,iii))=[];
                        
                    end
                    
                end
                
            end
            
        end
        
    end

    function do_predetermined()
        
        if ~isempty(pred_vars)
            
            just_neat=@myreplace; %#ok<NASGU>
            
            % add {0} to guys that are current
            %----------------------------------
            xpress0=parser.cell2matize(pred_vars);
            xpress=['\<',xpress0,'\>(?!{)'];
            model_eqtns=regexprep(model_eqtns,xpress,'$1{0}');
            
            % substract 1 to the time indices
            %---------------------------------
            xpress=['\<',xpress0,'\>{(\+|\-)?(\d+)}'];
            repl='${just_neat($1,$2,$3)}';
            model_eqtns=regexprep(model_eqtns,xpress,repl);
            
            % remove {0}
            %------------
            model_eqtns=strrep(model_eqtns,'{0}','');
            
        end
        
        function b=myreplace(a1,a2,a3)
            
            b=a1;
            
            a2a3_1=str2double([a2,a3])-1;
            
            if a2a3_1
                
                b=[b,'{',int2str(a2a3_1),'}'];
                
            end
            
        end
        
    end

    function add_stdev_to_model()
        
        expre=cell2mat(strcat(exo_list,'|'));
        
        expre=['\<(',expre(1:end-1),')\>'];
        
        repl=[stderr_name,'_$1*$1'];
        
        model_eqtns=regexprep(model_eqtns,expre,repl);
        
    end

    function rise_code = recreate_code()
        
        rise_code = [
            endo_block,eol,eol,...
            exo_block,eol,eol,...
            param_block,eol,eol];
        
        if ~isempty(obs_block)
            
            rise_code = [rise_code,...
                obs_block,eol,eol];
            
        end
        
        rise_code = [rise_code,...
            'model ',eol,eol,model_eqtns,eol,eol
            ];
        
        if ~isempty(ssmodel_eqtns)
            
            rise_code = [rise_code,...
                'steady_state_model ',eol,eol,ssmodel_eqtns,eol,eol
                ];
            
        end
        
        if ~isempty(planner_block)
            
            rise_code = [rise_code,eol,eol,...
                planner_block,eol,eol
                ];
        end
        
        rise_code = strrep(rise_code,';',sprintf(';\n\n'));
        
        rise_code = [
            regexprep(['% ',header],'(\n)','$1%'),eol,eol,...
            rise_code
            ];
        
        % remove the semicolon in lines starting with @#
        %-----------------------------------------------------------
        rise_code = regexprep(rise_code,'(@\s*#[^\n]+);','$1');
        
        
        parser.write2file(rise_code,riseFileName);
        
    end

    function list=extract_other_blocks(typeof,do_replace_time)
        
        if nargin<2
            
            do_replace_time=false;
            
        end
        
        tokens = regexpi(rise_code,['\<',typeof,'\s*(.*?)end;'],'tokens');
        
        if isempty(tokens)
            
            list=tokens;
            
            return
            
        end
        
        tokens=[tokens{:}];
        
        for iblk=1:numel(tokens)-1
            
            tokens{iblk}=[tokens{iblk},' '];
            
        end
                
        list = cell2mat(tokens);
        
        % replace time indices
        if do_replace_time
            
            list = regexprep(list,'(?<=\w)\(([\+-]?\d*)\)','{$1}');
            
        end
        
        list=regexprep(list,'(;|¤)\s+','$1');
        
        % add empty line at end of @#...
        patt='(@#\s*[^¤]+\s*)¤'; repl='$1\n\n';
        
        list=regexprep(list,patt,repl);
        
        list = strrep(list,';¤',';');
        
    end

    function insert_all_subfiles()
        
        patt='@#\s*include\s*"([^"]+)"';
        
        while true
            
            mm=regexp(raw_code,patt,'match');
            
            if isempty(mm)
                
                break
                
            end
            
            if ischar(mm)
                
                mm={mm};
                
            end
            
            for ifile=1:numel(mm)
                
                process_import(mm{ifile});
                
            end
            
        end
        
        function process_import(include_file)
            
            dblq=strfind(include_file,'"');
            
            fname=include_file(dblq(1)+1:dblq(2)-1);
            
            c= read_file(fname);
            
            lc=length(c);
            
            lr=length(raw_code);
            
            raw_code=strrep(raw_code,include_file,[eol,c,eol]);
            
            fprintf(1,'Uploading %s : old(%0.0f)+new(%0.0f)= %0.0f\n',...
                include_file,lr,lc,lr+lc);
            
        end
        
    end

end

function raw_code = read_file(dynFileName)

fid = fopen(dynFileName,'r');

if fid == -1
    
    if exist(dynFileName,'file') == false
        
        error('Unable to find ''%s''.',dynFileName);
        
    else
        
        error('Unable to open ''%s'' for reading.',dynFileName);
        
    end
    
end

raw_code = transpose(fread(fid,'*char'));

fclose(fid);

end

function est=extract_estimation(estim_block,std_type,obs_list)

% remove rows starting with corr
corr_rows=regexp(estim_block,'corr[^;]+;','match');

estim_block=regexprep(estim_block,'corr[^;]+;','');

% remove the trigger
estim_block=strrep(estim_block,'estimated_params;','');
estim_block=strrep(estim_block,'end;','');

splits=regexp(estim_block,';','split');

good=cellfun(@(x)~isempty(x),splits,'uniformOutput',true);

splits=splits(good);

n=numel(splits);

est=struct('name',{},'start',{},'lb',{},'ub',{},'distr',{},...
    'mean',{},'sd',{},'p3',{},'p4',{},'scale',{},'is_idistr_shift',{});

fields=fieldnames(est);

est0=cell(n,10);

% shiftBoundsDistr={'normal','gamma','inv_gamma','uniform','beta'};


for irow=1:n
    
    rawline=regexp(splits{irow},',','split');
    
    rawline=cellfun(@(x)x(~isspace(x)),rawline,'uniformOutput',false);
    
    est0(irow,1)=rawline(1);
    
    rawline(1)=[];
    
    iter=1;
    
    if ~isempty(rawline{1}) && isstrprop(rawline{1}(1),'alpha')
        
        iter=5;
        
    else
        
        iter=iter+1;
        
    end
    
    while ~isempty(rawline)
        
        est0{irow,iter}=rawline{1};
        
        rawline(1)=[];
        
        iter=iter+1;
        
    end
    
    if isempty(est0{irow,2})
        
        est0{irow,2}=est0{irow,6};
        
    end
    
    if isempty(est0{irow,5})
        
        est0{irow,5}='uniform';
        
    else
        
        est0{irow,5}=strrep(lower(est0{irow,5}),'_pdf','');
        
    end
    
    est(irow).name=est0{irow,1};
    
    est(irow).start=est0{irow,2};
    
    est(irow).lb=est0{irow,3};
    
    est(irow).ub=est0{irow,4};
    
    est(irow).distr=est0{irow,5};
    
    est(irow).mean=est0{irow,6};
    
    est(irow).sd=est0{irow,7};
    
    est(irow).p3=est0{irow,8};
    
    est(irow).p4=est0{irow,9};
    
    est(irow).scale=est0{irow,10};
    
    est(irow).is_idistr_shift=false;
    
    est(irow)=process_one(est(irow));
    
end

shifted=est([est.is_idistr_shift]);

if ~isempty(shifted)
    
    warning(['RISE does not support beta distributions ',...
        'with ranges outside [0,1]. Switching to a truncated ',...
        'normal distribution for parameter(s): '])
    
    disp({shifted.name})
    
end


    function x=process_one(x)
        
        x.name=process_name(x.name);
        
        for ifield=1:numel(fields)
            
            ff=fields{ifield};
            
            if ~isempty(x.(ff))
                
                x.(ff)(isspace(x.(ff)))=[];
                
            end
            
        end
        
        if ~isempty(x.p3) % any(strcmp(x.distr,shiftBoundsDistr)) && 
            
            x.lb=x.p3;
            
            if ~isempty(x.p4)
                
                x.ub=x.p4;
                
            end
            
        end
        
        lb=str2double(x.lb);
        
        ub=str2double(x.ub);
        
        if ~isnan(lb)
            
            if (strcmp(x.distr,'gamma') && lb==0 && ub==inf)||...
                    (strcmp(x.distr,'beta') && lb==0 && ub==1)
                
                x.lb=[];
                
                x.ub=[];
                
                
            elseif (strcmp(x.distr,'beta') && ~((lb>=0 && lb<=1) && (ub<=1)))
                
                x.is_idistr_shift=true;
                
                x.distr='normal';
                
            end
            
            
        end
        
        if (strcmp(x.distr,'inv_gamma') && strcmpi(x.sd,'inf'))
            
            warning(['RISE does not understand inf ',...
                'changing the prior standard deviation of parameter ',...
                x.name,' from inf to 4'])
            
            x.sd='4';
            
        end
        
        function name=process_name(name)
            % get rid of any space before the name
            name=regexp(name,'\w+.*','match');
            
            name=name{1};
            
            if strncmp(name,'stderr',6)
                
                name=strrep(name,'stderr','');
                
                name(isspace(name))=[];
                
                main_type=std_type;
                
                if ismember(name,obs_list)
                    
                    main_type='stderr';
                    
                end
                
                name=[main_type,'_',name];
                
            end
            
        end
        
    end

end

function rise_code=replace_descriptions(rise_code)
% replace y ${y}$ (long_name='output') with y "{y}(output)"
%----------------------------------------------------------
patt2='(\w+[^$]*)\s*,?\s*\$(.*?)\$\s*,?\s*\(\s*long_name\s*=\s*''([^'']+)''\s*\)';
% patt2='(\w+)\s*,?\s*\$(.*?)\$\s*,?\s*\(\s*long_name\s*=\s*''([^'']+)''\s*\)';
% replace2='$1 "$3($2)"';
replace2='$1 "$3 # $2"';

replacer=@replace_engine; %#ok<NASGU>
patt='\<(var|varexo|parameters)\>([^;]+;)';
repl='${replacer($1,$2)}';
rise_code = regexprep(rise_code,patt,repl);

    function out=replace_engine(str1,str2)
        
        out = regexprep(str2,patt2,replace2);
        
        out=[str1,out];
        
    end

end