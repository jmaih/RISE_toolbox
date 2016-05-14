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
raw_code = read_file();

% remove block comments
%-----------------------
rise_code = regexprep(raw_code,'/\*(.*)\*/','');

header = sprintf('Conversion of Dynare file [%s] into RISE file [%s]\n',...
    dynFileName,riseFileName);

% Remove line comments.
%----------------------
rise_code = regexprep(rise_code,'(//|%)(.*?\n)','\n');

% replace y ${y}$ (long_name='output') with y "{y}(output)"
%----------------------------------------------------------
express='(\w+)\s*,?\s*\$(.*?)\$\s*,?\s*\(\s*long_name\s*=\s*''([^'']+)''\s*\)';
% replace='$1 "$3($2)"';
replace='$1 "$3 # $2"';
rise_code = regexprep(rise_code,express,replace);

% % add a semicolon at the end of each line starting with a @#, it will be
% % removed later
% %--------------------------------------------------------------------------
% rise_code = regexprep(rise_code,'(@\s*#[^\n;]+)','$1;');

% place a ¤ at the end of each @#...
%-------------------------------------------
atPound='@\s*#\s*[^\n]+';
express=['(',atPound,')\n'];
rise_code = regexprep(rise_code,express,'¤$1¤');

% replace all $ with "
%---------------------
rise_code=strrep(rise_code,'$','"');

% replace all endif with end
%----------------------------
rise_code=regexprep(rise_code,'(#\s*)\<endif\>','$1end');

% add space to all # signs
%----------------------------
rise_code=strrep(rise_code,'#','# ');

% Convert char(10) to white space
%--------------------------------
eol=char(10);% sprintf('\n')
tab=char(9);% sprintf('\t')

% remove multiple white characters
%----------------------------------
rise_code = regexprep(rise_code,'\s+',' ');

% predetermined variables
%-------------------------
[~,pred_vars]=extract_declaration_block('predetermined_variables');

% extract block of endogenous
%-----------------------------
[endo_block]=extract_declaration_block('var','endogenous');

% extract block of exogenous
%-----------------------------
[exo_block,exo_list]=extract_declaration_block('varexo','exogenous');

% extract block of exogenous
%-----------------------------
[param_block,par_list]=extract_declaration_block('parameters');

% extract block of observables
%-----------------------------
[obs_block]=extract_declaration_block('varobs','observables');

% extract planner objective
%---------------------------
[planner_block]=extract_planner_objective();

% add new parameters at the end of the parameters block. Those parameters
% will not have definitions. But what if the block is empty?
%-----------------------------------------------------------------------
newparams=cell2mat(strcat([stderr_name,'_'],exo_list,'@'));

newparams=strrep(newparams,'@',' ');

param_block=[param_block,' ',newparams];

% model equations
%-----------------
rise_code = regexprep(rise_code,'model\(linear\)','model');
model_eqtns = extract_other_blocks('model;',true);

% take care of predetermined variables
do_predetermined();

% steady state model equations
%------------------------------
ssmodel_eqtns = extract_other_blocks('steady_state_model;');

% push standard deviations into equations directly
%-------------------------------------------------
add_stdev_to_model()

% Looking for covariances: shocks section
%----------------------------------------
shocks_block = extract_other_blocks('shocks;');

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
        
        str=regexp(rise_code,express,'names');
                
        pnames={str.name};
        
        % taking care of recursive computations
        xpress=['\<',parser.cell2matize(pnames),'\>'];
        
        repl='p.$1';
        
        pvals=regexprep({str.value},xpress,repl);
        
        % add shock standard deviations
        shock_standard_deviations();
        
        pvals=cellfun(@(x)x(~isspace(x)),pvals,'uniformOutput',false);
        
        code=strcat('p.',pnames,'=',pvals,';',sprintf('\n'));
        
        paramFileName=regexprep(riseFileName,'(\w+)\.\w+','$1_params');
        
        manualAdjust=['% Remarks: ',sprintf('\n'),...
            '% - Only the parameters and shocks variances or ',sprintf('\n'),...
            '%   standard deviations found in the dynare file are assigned. ',...
            sprintf('\n'),...
            '% - RISE will set all parameters without a value to nan'];
        
        code=[sprintf('function p=%s()\n%s\n\np=struct();\n',...
            paramFileName,manualAdjust),...
            code];
        
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

    function [blk,list]=extract_declaration_block(trigger,repl_trigger)
        
        express_=['(\<',trigger,'\>\s+[^;]+;)'];
        
        blk=regexp(rise_code,express_,'tokens','once');
        
        if nargin>1
            
            blk=regexprep(blk,['\<',trigger,'\>'],repl_trigger);
            
        else
            
            repl_trigger=trigger;
            
        end
        
        % remove the semicolon
        %----------------------
        if isempty(blk)
            
            blk='';
            
            list={};
            
            return
            
        end
        
        blk=blk{1};
        
        semcol=find(blk==';',1,'last');
        
        blk(semcol)=[];
        
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
            
            list=regexp(tmp,'\w+','match');
            
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
            regexprep(['%',header],'(\n)','$1%'),eol,eol,...
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
        
        tokens = regexpi(rise_code,[typeof,'\s*(.*?)end;'],...
            'tokens','once');
        
        if isempty(tokens)
            
            list=tokens;
            
            return
            
        end
        
        list = tokens{1};
        
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

    function raw_code = read_file()
        
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

end

