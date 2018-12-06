function duplicates=dynare2rise(dynFileName,riseFileName,...
    stderr_name,detrend)
% dynare2rise -- converts a basic dynare file into a RISE one
%
% ::
%
%   duplicates=dynare2rise(dynFileName)
%   duplicates=dynare2rise(dynFileName,riseFileName)
%   duplicates=dynare2rise(dynFileName,riseFileName,stderr_name)
%   duplicates=dynare2rise(dynFileName,riseFileName,stderr_name,detrend)
%
% Args:
%
%    - **dynFileName** [char] : name of the dynare model file with or
%      without extension (.dyn or .mod)
%
%    - **riseFileName** [char|{'dynFileName.rs'}] : name of the created
%      RISE file 
%
%    - **stderr_name** [char|{'std'}] : prepended name for the newly
%      created parameters (RISE transforms all the variances and standard
%      deviations into parameters)
%
%    - **detrend** [true|{false}]: 
%      - when true, the model is written in stationary form like dynare
%        would do behind the scenes. For log growth variables,
%        y-->(y*trend). For linear growth variables, y-->(y+trend). The
%        steady state model in dynare then corresponds exactly to the
%        steady state file produced by RISE, except that the trend
%        variables are added as new variables. 
%      - when false, the model is kept in nonstationary form. For log
%        growth variables the growth rate in an expression of the form f(x)
%        becomes log(f(exp(x)). The growth rate for linear growth variables
%        remains unchanged. The variables detected to have steady state = 0
%        but that need to be stationarized are assigned a growth rate of 0
%        and are declared as LEVEL variables. The steady state remains
%        unchanged.
%      - The steady state for all trend variables with log growth is
%        normalized to 1, while the steady state for trend variables with
%        linear growth is normalized to 0.
%
% Returns:
%    :
%
%    - **duplicates** [struct] : structure containing information on the
%      redundancies
%
% NB: the created files are 
%
%    - riseFileName.rs 
%
%    - riseFileName_params.m : parameter function file with calibration and
%      priors as separate outputs
%
%    - riseFileName_sstate.m (optional): steady state function file
%
% See also :


if nargin<4
    
    detrend=[];
    
    if nargin<3
        
        stderr_name=[];
        
        if nargin<2
            
            riseFileName=[];
            
        end
        
    end
    
end

duplicates=struct();

if isempty(detrend)
    
    detrend=false;
    
end

paramFileName='';

if isempty(riseFileName)
    
    riseFileName=[regexprep(dynFileName,'(\w+)\.\w+','$1'),'.rs'];
    
else
    
    xt=find(riseFileName=='.',1);
    
    if isempty(xt)
        
        riseFileName=[riseFileName,'.rs'];
        
    end
    
end

if isempty(stderr_name)
    
    stderr_name='std';
    
end

dynFileName = strtrim(dynFileName);

riseFileName = strtrim(riseFileName);

if exist('contains.m','file')
    % a contains b
    my_contains=@(a,b)contains(a,b);
    
else
    % find b in a
    my_contains=@(a,b)~isempty(strfind(a,b)); %#ok<STREMP>
    
end

% read input file
%-----------------
raw_code = read_file(dynFileName,my_contains);

% Convert char(10) to white space
%--------------------------------
if exist('newline.m','file')
    
    eol=newline();
    
else
    
    eol=char(10);%#ok<CHARTEN> % sprintf('\n')
    
end
% tab=char(9);% sprintf('\t')

% insert all subfiles
%--------------------
insert_all_subfiles();

header = sprintf('Conversion of Dynare file [%s] into RISE file [%s]\n',...
    dynFileName,riseFileName);

rise_code=raw_code;

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

% extract trend vars
%-------------------
trend_var=struct();
log_trend_var=struct();
extract_trend_var('trend_var');
extract_trend_var('log_trend_var');

% extract block of endogenous
%-----------------------------
is_endo_decl=true;
deflators=struct();
log_deflators=struct();
[endo_block,endo_list]=extract_declaration_block('var','endogenous',is_endo_decl);

% add the trend_var to the list
%-------------------------------
trend_list=fieldnames(trend_var);
log_trend_list=fieldnames(log_trend_var);
endo_list=union(endo_list,trend_list);
endo_list=union(endo_list,log_trend_list);

% collect log variables
%-----------------------
log_variables=trend_list;

level_variables=log_trend_list;

if ~detrend
    
    log_variables=union(fieldnames(deflators),log_variables);
    
    level_variables=union(fieldnames(log_deflators),log_trend_list);
    
end

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

% new parameters in a separate list that will be used later on
%-------------------------------------------------------------
new_par_list=regexp(newparams,'\w+','match');

% complete the list of parameters
%--------------------------------
par_list=[par_list,new_par_list];

% model equations: discard attributes
%--------------------------------------
rise_code = regexprep(rise_code,'model\([^)]+\)','model');

rise_code = regexprep(rise_code,'(\+\s*\-|\-\s*\+)','-');

do_replace_time=false;
model_eqtns = extract_other_blocks('model',do_replace_time);

do_model_conversions()

% steady state model equations
%------------------------------
ssmodel_eqtns = extract_other_blocks('(steady_state_model|initval)');

if ~isempty(ssmodel_eqtns)
    
    ssmodel_eqtns=strrep(ssmodel_eqtns,'initval','');
    
    ssmodel_eqtns=strrep(ssmodel_eqtns,'steady_state_model','');
    
end

% balanced growth path model
%---------------------------
bgp_model_eqtns='';

% add the steady state for the trend vars to the steady state
% equations
%-------------------------------------------------------------
ss_trend_var=sstate_for_trend_variables(trend_var,'trend');

ss_log_trend_var=sstate_for_trend_variables(log_trend_var,'log_trend');

process_group_of_deflators(deflators,'trend')

process_group_of_deflators(log_deflators,'log_trend')

bgp_model_eqtns(isspace(bgp_model_eqtns))=[];

ssmodel_eqtns=[
    '    % trend variables ',eol,...
    '    %------------------',eol,...
    ss_trend_var,...
    '    % log trend variables ',eol,...
    '    %---------------------',eol,...
    ss_log_trend_var,...
    '    % old steady state ',eol,...
    '    %------------------',eol,...
    ssmodel_eqtns];

% Looking for covariances: shocks section
%----------------------------------------
shocks_block = extract_other_blocks('shocks');

% Looking for covariances: shocks section
%----------------------------------------
estim_block = extract_other_blocks('estimated_params');

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

% list of parameters assigned during steady state calculation
newp_list=[];
% list of endogenous variables solved in steady state
ssList=[];

write_parameter_file()

write_steady_state_file()

recreate_code();

untagged_model=untag_model(model_eqtns);

neqtns=sum(untagged_model==';');

summary=[' ',eol,...
    ['# equations : ',int2str(neqtns)],eol,...
    ['# endogenous : ',int2str(numel(endo_list))],eol,...
    ['# endogenous solved in sstate: ',int2str(numel(ssList))],eol,...
    ['# endogenous(duplicates) solved in sstate: ',int2str(numel(duplicates.ssList))],eol,...
    ['# exogenous : ',int2str(numel(exo_list))],eol,...
    ['# parameters : ',int2str(numel(par_list))],eol,...
    ['# parameters assigned in sstate: ',int2str(numel(newp_list))],eol,...
    ['# parameters(duplicates) assigned in sstate: ',int2str(numel(duplicates.newp_list))],eol,...
    ['# parameters created for RISE: ',int2str(numel(new_par_list))],eol,...
    ['# parameters(original) : ',int2str(numel(par_list)-numel(new_par_list))]
    ];

fprintf('%s\n',summary);

    function process_group_of_deflators(batch,type_)
        
        cylist=parser.cell2matize(endo_list(:).');
        
        list=fieldnames(batch);
        
        switch type_
            
            case 'trend'
                
                is_linear_trend=false;
                
            case 'log_trend'
                
                is_linear_trend=true;
                
            otherwise
                
                error('This should never happen')
                
        end
        
        for iii=1:numel(list)
            
            vv=list{iii};
            
            vg=batch.(vv);
            
            if ~detrend
                
                bgp_model_eqtns=[bgp_model_eqtns,...
                    vv,'=',process_one_deflator(vg),';']; %#ok<AGROW>
                
            end
            
        end
        
        function vg=process_one_deflator(vg0)
            
            if all(isstrprop(vg0,'alphanum'))||is_linear_trend
                
                vg=vg0;
                
            else
                
                vg=regexprep(vg0,['\<',cylist,'\>'],'exp($1)');
                
                vg=['log(',vg,')'];
                
            end
            
        end
        
    end

    function sstates=sstate_for_trend_variables(batch,type_)
        
        list=fieldnames(batch);
        
        if strcmp(type_,'trend')
            
          engine=@(x)['log(',batch.(x),')'];
          
          ss='1';
          
        elseif strcmp(type_,'log_trend')
            
          engine=@(x)batch.(x);
          
          ss='0';
          
        else
            
            error('This should not happen')
            
        end
        
        sstates='';
        
        for iii=1:numel(list)
            
            vv=list{iii};
            
            sstates=[sstates,'    ',vv,'=',ss,';']; %#ok<AGROW>
            
            bgp_model_eqtns=[bgp_model_eqtns,'    g.',...
                vv,'=',engine(vv),';']; %#ok<AGROW>
            
        end
        
    end

    function write_steady_state_file()
        
        % make list of new parameters
        tmp=parser.cell2matize(par_list(:).');
        
        newp_list=regexp(ssmodel_eqtns,['\<',tmp,'\>\s*='],'tokens');
        
        [newp_list,duplicates.newp_list]=set_list([newp_list{:}]);
        
        olp_list=par_list-newp_list;
        
        yList=parser.cell2matize(endo_list(:).');
        
        olp_list__=parser.cell2matize(olp_list(:).');
        
        newp_list__=parser.cell2matize(newp_list(:).');
        
        write_ss_model_eqtns=regexprep(ssmodel_eqtns,['\<',yList,'\>'],'v.$1');
        
        % locate the variables that are pertinent
        %-----------------------------------------
        ssList=regexp(write_ss_model_eqtns,'\<v\.(\w+)\s*=','tokens');
        
        [ssList,duplicates.ssList]=set_list([ssList{:}]);
        
        % Declared variables with exactly zero steady state
        %---------------------------------------------------
        declared_ss0_vars_List=regexp(write_ss_model_eqtns,'\<v\.(\w+)\s*=\s*0\s*;','tokens');
        
        [declared_ss0_vars_List,duplicates.declared_ss0_vars_List]=set_list([declared_ss0_vars_List{:}]);
        
        non_ss0_list=ssList-declared_ss0_vars_List;
        
        % all variables with zero steady state
        %--------------------------------------
        all_ss0_vars_List=endo_list-non_ss0_list;
        
        % balanced-growth path equations
        %--------------------------------                
        write_ss_model_eqtns=regexprep(write_ss_model_eqtns,['\<',newp_list__,'\>'],'newp.$1');
        
        write_ss_model_eqtns=regexprep(write_ss_model_eqtns,['\<',olp_list__,'\>'],'p.$1');
        
        write_bgp_model_eqtns=regexprep(bgp_model_eqtns,['(?<!\.)\<',yList,'\>'],'g.$1');
        
        write_bgp_model_eqtns=regexprep(write_bgp_model_eqtns,['\<',newp_list__,'\>'],'newp.$1');
        
        write_bgp_model_eqtns=regexprep(write_bgp_model_eqtns,['\<',olp_list__,'\>'],'p.$1');
        
        if ~detrend
            % ss = 0 variables with non-zero bgp
            %-----------------------------------
            tmp=parser.cell2matize(all_ss0_vars_List(:).');
            
            ss0_bgp_list=regexp(write_bgp_model_eqtns,['g.',tmp,'s*=[^;]+;'],'tokens');
            
            [ss0_bgp_list,duplicates.ss0_bgp_list]=set_list([ss0_bgp_list{:}]);
            
            % remove (Comment out) all bgp variables that have sstate = 0
            %-------------------------------------------------------------
            % write_bgp_model_eqtns=regexprep(write_bgp_model_eqtns,['g\.',tmp,'s*=[^;]+;'],'');
            write_bgp_model_eqtns=regexprep(write_bgp_model_eqtns,['(g\.)',tmp],'%$1$2');
            
            log_variables=log_variables-ss0_bgp_list;
            
        end
        
        code=[
            'retcode=0;',eol,eol,...
            ['ylist = ',make_list(ssList),';'],eol,eol,...
            '% ylist = get(m,''endo_list(original)'');',eol,eol,...
            'if nargin==1',eol,...
            '    % list of endogenous variables to be calculated',eol,...
            '    %----------------------------------------------',eol,...
            '    y=ylist;',eol,...
            '    % list of parameters to be computed during steady state calculation',eol,...
            '    %-------------------------------------------------------------------',eol,...
            ['    newp=',make_list(newp_list),';'],eol,...
            'else',eol,...
            '    % steady states and new parameter values',eol,...
            '    %---------------------------------------',eol,...
            '    v=struct();',eol,eol,...
            '    newp = struct();',eol,eol,...
            '    % dynare implicitly assumes all variables to take 0 at steady state',eol,eol,...
            '    y(:,1) = 0;',eol,eol,...
            strrep(write_ss_model_eqtns,';',[';',eol,eol]),eol,eol,...
            '    % growth rates',eol,...
            '    %-------------',eol,...
            '    g = struct();',eol,...
            strrep(write_bgp_model_eqtns,';',[';',eol,eol]),eol,eol,...
            '    N = length(ylist);',eol,eol,...
            '    is_nonstationary=size(y,2)==2;',eol,eol,...
            '    ys=zeros(length(ylist),1+is_nonstationary);',eol,eol,...
            '    for jn=1:N',eol,eol,...
            '       if isfield(v,ylist{jn})',eol,eol,...
            '           ys(jn,1) = v.(ylist{jn});',eol,eol,...
            '       end',eol,eol,...
            '       if is_nonstationary && isfield(g,ylist{jn})',eol,eol,...
            '           ys(jn,2) = exp(g.(ylist{jn}));',eol,eol,...
            '       end',eol,eol,...
            '    end',eol,eol,...
            '',eol,...
            '    % check the validity of the calculations',eol,...
            '    %----------------------------------------',eol,...
            '    if ~utils.error.valid(ys)',eol,eol,...
            '        retcode=1;',eol,eol,...
            '    else',eol,...
            '',eol,...
            '        % push the calculations',eol,...
            '        %----------------------',eol,...
            '        y(id,:)=ys;',eol,eol,...
            '    end',eol,eol,...
            'end',eol,eol
            ];
        
        ssfile_name=regexprep(riseFileName,'(\w+)\.\w+','$1_sstate');
        
        manualAdjust=[
            ['% ',ssfile_name,'--  computes the steady state of the model analytically '],eol,...
            '% ',eol,...
            '% ::',eol,...
            '% ',eol,...
            '% ',eol,...
            ['%   [y,newp,retcode]=',ssfile_name,'(m,y,p,d,id)'],eol,...
            '% ',eol,...
            '% Args:',eol,...
            '% ',eol,...
            '%   - **m** [rise|dsge]: model object (not always needed)',eol,...
            '% ',eol,...
            '%   - **y** [vector]: endo_nbr x 1 vector of initial steady state ',eol,...
            '% ',eol,...
            '%    - **pp** [struct]: parameter structure ',eol,...
            '% ',eol,...
            '%    - **d** [struct]: definitions ',eol,...
            '% ',eol,...
            '%    - **id** [vector]: location of the variables to calculate ',eol,...
            '% ',eol,...
            '% Returns: ',eol,...
            '%     :',eol,...
            '% ',eol,...
            '%    - **y** []: endo_nbr x 1 vector of updated steady state ',eol,...
            '% ',eol,...
            '%    - **newp** [struct]: structure containing updated parameters if any ',eol,...
            '% ',eol,...
            '%    - **retcode** [0|number]: return 0 if there are no problems, else return ',eol,...
            '%      any number different from 0 ',eol,...
            '% ',eol,...
            '% Note: ',eol,...
            '% ',eol,...
            '%    - this is new approach has three main advantages : ',eol,...
            '%      - The file is valid whether we have many regimes or not ',eol,...
            '%      - The user does not need to know what regime is being computed ',eol,...
            '%      - It is in sync with the steady state model ',eol,...
            '% ',eol,...
            '% Example: ',eol,...
            '% ',eol,...
            '%    See also: ',eol,...
            '% ',eol,...
            regexprep(['% ',header],'(\n)','$1%'),eol,...
            '% ',eol];
        
        code=[sprintf('function [y,newp,retcode]=%s(m,y,p,d,id) %%#ok<INUSL>\n%s\n\n',...
            ssfile_name,manualAdjust),...
            code,eol,'end'];
        
        parser.write2file(code,[ssfile_name,'.m'])
        
        function list=make_list(list)
            
            if isempty(list)
                
                list='{}';
                
            else
                
                list=cell2mat(strcat(list(:).',','));
                
                list=['''',strrep(list(1:end-1),',',''','''),''''];
                
                list=['{',list,'}'];
                
            end
            
        end
        
    end

    function write_parameter_file()
        
        % parameter values
        %------------------
                
        % convert set_param_value('pname',expression); into
        % pname=expression;
        rise_code=regexprep(rise_code,'\<set_param_value\s*\(\s*''(\w+)''\s*,([^)]+)\);','$1=$2;');
        
        % Do the following sweep at least twice in order to collect
        % temporary/local/auxiliary names that are not defined as
        % parameters but that are used to define parameters
        auxil_params={};
        
        while 1
            % Do not parse parameters that are preceded with a "." they are
            % most likely matlab...
            express_=['(?<!\.)(?<name>\<',parser.cell2matize(par_list),'\>)\s*=(?<value>[^;]+);'];
            
            str=regexp(rise_code,express_,'names');
            
            pnames={str.name};
            
            pvals={str.value};
            
            % Add all atoms that are not parameters or matlab functions to the
            % list of parameters
            tmp=regexp(pvals,'\<[a-zA-Z]+\w*\>(?!\()','match');
            
            tmp=[tmp{:}];
            
            tmp=tmp-par_list;
            
            if isempty(tmp)
                
                break
                
            end
            
            disp(tmp)
            
            par_list=[par_list,tmp]; %#ok<AGROW>
            
            auxil_params=tmp;
            
        end
        
%         express_=['\<',parser.cell2matize(tmp),'\>'];
%         pvals=regexprep(pvals,express_,'p.$1');
        
        % add shock standard deviations
        shock_standard_deviations();
        
        % taking care of recursive computations: use all parameters and not
        % just those that were assigned in the model file
        xpress=['\<',parser.cell2matize(par_list),'\>']; % <-- xpress=['\<',parser.cell2matize(pnames),'\>'];

        repl='p.$1';
        
        pvals=regexprep(pvals,xpress,repl);
        
        pvals=cellfun(@(x)x(~isspace(x)),pvals,'uniformOutput',false);
        
        myditch=@ditch_redundant; %#ok<NASGU>
        
        code=strcat('p.',pnames,'=',pvals,';',eol);
        
        code=regexprep(code,'(p.\w+)\s*=\s*([^;]+);','${myditch($1,$2)}');
        % Now undo the move above: remove auxiliary parameters
        %-----------------------------------------------------
        if ~isempty(auxil_params)
            
            par_list=par_list-auxil_params;
            
            xp=['p.\<',parser.cell2matize(auxil_params),'\>'];
            
            code=regexprep(code,xp,'$1');
            
        end
                               
        code=['if nargin>1',eol,...
            '   p=load(matfile);',eol,...
            'else',eol,...
            '   p=struct();',eol,...
            'end',eol,...
            code,eol,...
            '% discarding atoms not declared as parameters in the model',...
            '%---------------------------------------------------------',...
            'np=struct();',eol,...
            'if nargin && ~isempty(m) && isa(m,''dsge'')',eol,...
            '   fp=fieldnames(p);',eol,...
            '   plist=get(m,''par_list'');',eol,...
            '   bad=~ismember(fp,plist);',eol,...
            '   tmp=rmfield(p,fp(bad));',eol,...
            '   np=rmfield(p,fp(~bad));',eol,...
            '   p=tmp;',eol,...
            'end'
            ];
        
        paramFileName=regexprep(riseFileName,'(\w+)\.\w+','$1_params');
        
        manualAdjust=[['% ',paramFileName,'--  Sets the baseline calibration and the priors for estimation '],eol,...
            '% ',eol,...
            '% ::',eol,...
            '% ',eol,...
            '% ',eol,...
            ['%   [p,priors,np]=',paramFileName,'(m)'],eol,...
            ['%   [p,priors,np]=',paramFileName,'(m,matfile)'],eol,...
            '% ',eol,...
            '% Args:',eol,...
            '% ',eol,...
            '%   - **m** [rise|dsge]: model object (optional)',eol,...
            '%     only used to discard parameters that are not declared in the model',eol,...
            '%   - **matfile** [.mat|{[]}]: .mat file containing some other parameterization  (optional)',eol,...
            '%     only used to discard parameters that are not declared in the model',eol,...
            '% ',eol,...
            '% Returns: ',eol,...
            '%     :',eol,...
            '% ',eol,...
            '%    - **p** [struct]: structure containing parameter values ',eol,...
            '% ',eol,...
            '%    - **priors** [struct]: structure containing the priors for estimation in RISE ',eol,...
            '% ',eol,...
            '%    - **np** [struct]: structure containing atoms not declared as parameters in the model ',eol,...
            '% ',eol,...
            '% Note: ',eol,...
            '% ',eol,...
            '% - The parameters and shocks variances or ',eol,...
            '%   standard deviations found in the dynare file are assigned. ',eol,...
            '% - Shock standard deviations not assigned in the ',eol,...
            '%   dynare file get a value of 0 following dynare''s convention',eol,...
            '% - RISE will set all other parameters without a value to nan',eol,...
            '% ',eol,...
            '% Example: ',eol,...
            '% ',eol,...
            '%    See also: ',eol,...
            '% ',eol,...
            regexprep(['% ',header],'(\n)','$1%'),eol,...
            '% ',eol
            ];
        
        codest=[
            '% priors for estimation',eol,...
            '%----------------------',eol,...
            'priors=struct();',eol,eol
            ];
        
        for ip=1:numel(est)
            
            par_i=est(ip);
            
            if strcmp(par_i.distr,'uniform')
                
                my_new_line=sprintf('priors.%s={%s, %s, %s',...
                    par_i.name,par_i.start,par_i.lb,par_i.ub);
                
            else
                
                my_new_line=sprintf('priors.%s={%s, %s, %s, ''%s''',...
                    par_i.name,par_i.start,par_i.mean,par_i.sd,par_i.distr);
                
                if ~isempty(par_i.lb)
                    
                    my_new_line=sprintf('%s, %s',my_new_line,par_i.lb);
                    
                    if ~isempty(par_i.ub)
                        
                        my_new_line=sprintf('%s, %s',my_new_line,par_i.ub);
                        
                    end
                    
                end
                
            end
            
            codest=[codest,my_new_line,'};',eol,eol]; %#ok<AGROW>
            
        end
        
        code=[sprintf('function [p,priors,np]=%s(m,matfile)\n%s\n\n',...
            paramFileName,manualAdjust),...
            code,eol,codest];
        
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
        
        function out=ditch_redundant(left,right)
                
                out=[left,'=',right,';']; 
            
            if strcmp(left,right)
                
                out=['% ',out];
                
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
        
        isCommitment=my_contains(rise_code,'ramsey_policy');
        
        isDiscretion=false;
        
        if ~isCommitment
            
            isDiscretion=my_contains(rise_code,'discretionary_policy');
            
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
                
        if is_endo_decl
            % do not read from the var ... that are hidden inside the
            % shocks block
            % pull_from=regexprep(rise_code,'\<shocks\>;\s*(.*?)end;','');
            % alternatively, replace var with var__ inside the shocks
            % blocks and then undo the change after extraction
            blank_name='var_____';
            
            blank_var_in_shocks_blocks()
            
        end
        
        blk=regexp(rise_code,express_,'tokens');
        
        rise_code=regexprep(rise_code,express_,'');
        
        if is_endo_decl
            
            rise_code=strrep(rise_code,blank_name,'var');
            
        end
        
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
        
        bing_defl=regexp(blk,'(\s*deflator\s*=','start');
        
        bing_log_defl=regexp(blk,'(\s*log_deflator\s*=','start');
        
        for iblk=1:numel(blk)
            
            semcol=find(blk{iblk}==';',1,'last');
            
            blk{iblk}(semcol)=[];
            
            blk{iblk}=set_deflator(blk{iblk},bing_defl{iblk},'deflator');
            
            blk{iblk}=set_deflator(blk{iblk},bing_log_defl{iblk},'log_deflator');
            
            if iblk>1
                
                blk{iblk}=strrep(blk{iblk},repl_trigger,'');
                
            end
            
        end
        
        blk=cell2mat(blk);% blk=blk{1};
        
        if nargout>1
            
            list=extract_list_of_names();
            
        end
        
        blk=remove_declaration_anchors(blk);
        
        function item=set_deflator(item,start,type_)
            
            if isempty(start)
                
                return
                
            end
            
            nopen=true;
            
            for iii=start+1:length(item)
                
                if item(iii)=='('
                    
                    nopen=nopen+1;
                    
                elseif item(iii)==')'
                    
                    nopen=nopen-1;
                    
                end
                
                if nopen==0
                    
                    break
                    
                end
                
            end
            
            range=start:iii;
            
            stud=regexprep(item(range(1)+1:range(end)-1),[type_,'\s*=\s*'],'');
            
            item(range)=[];
            
            fake_item=item;
            
            dblq=find(fake_item=='"',2,'last');
            
            while ~isempty(dblq)
                
                fake_item(dblq(1):dblq(end))=[];
                
                dblq=find(fake_item=='"',2,'last');
                
            end
            
            tmp=regexp(fake_item,'\w+','match');
            
            is_deflator=strcmp(type_,'deflator');
            
            for iii=1:numel(tmp)
                
                if ~strcmp(tmp{iii},'endogenous')
                    
                    if is_deflator
                        
                        deflators.(tmp{iii})=stud;
                        
                    else
                        
                        log_deflators.(tmp{iii})=stud;
                        
                    end
                    
                end
                
            end
            
        end
        
        function x=remove_declaration_anchors(x)
            
            x=strrep(x,'¤',eol);
            
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
        
        function blank_var_in_shocks_blocks()
            
            locs=regexp(rise_code,'\<shocks\>\s*;\s*(.*?)end\s*;','tokenExtents');
            
            for iii=numel(locs):-1:1
                
                pos=locs{iii};
                
                batch=regexprep(rise_code(pos(1):pos(2)),'\<var\>',blank_name);
                
                rise_code=[rise_code(1:pos(1)-1),batch,rise_code(pos(2)+1:end)];
                
            end
            
        end
        
    end

    function rise_code = recreate_code()
        
        rise_code = [
            endo_block,eol,eol,...
            exo_block,eol,eol,...
            param_block,eol,eol];
        
        if ~isempty(log_variables)
            
            lv=cell2mat(strcat(log_variables(:).',','));
            
            log_vars=['log_variables ',lv(1:end-1)];
            
            rise_code = [rise_code,log_vars,eol,eol];
            
        end
        
        if ~isempty(obs_block)
            
            rise_code = [rise_code,...
                obs_block,eol,eol];
            
        end
        
        rise_code = [rise_code,...
            'model ',eol,eol,model_eqtns,eol,eol
            ];
                
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

    function extract_trend_var(type_)% trend_var
        
        growth_factor='growth_factor';
        
        is_log_growth=strncmp('log_',type_,4);
        
        if is_log_growth
            
            growth_factor=['log_',growth_factor];
            
        end
        
        expression=['\<',type_,'\>\s*(\s*',growth_factor,'\s*=\s*[^;]+;'];
        
        list=regexp(rise_code,expression,'match');
        
        for iii=1:numel(list)
            
            list{iii}=strrep(list{iii},';','');
            
            equal=find(list{iii}=='=');
            
            lp=find(list{iii}==')',1,'last');
            
            vname=list{iii}(lp+1:end);
            
            vname(isspace(vname))=[];
            
            if is_log_growth
                
                log_trend_var.(vname)=list{iii}(equal+1:lp-1);
                
            else
                
                trend_var.(vname)=list{iii}(equal+1:lp-1);
                
            end
            
        end
        % remove from code
        rise_code=regexprep(rise_code,expression,'');
        
    end

    function list=extract_other_blocks(typeof,do_replace_time)
        
        if nargin<2
            
            do_replace_time=false;
            
        end
        
        pattern=['\<',typeof,'\s*;\s*(.*?)end\s*;'];
        
        tokens = regexpi(rise_code,pattern,'tokens');
        
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
        
        rise_code = regexprep(rise_code,pattern,'');
        
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
            
            c= read_file(fname,my_contains);
            
            lc=length(c);
            
            lr=length(raw_code);
            
            raw_code=strrep(raw_code,include_file,[eol,c,eol]);
            
            fprintf(1,'Uploading %s : old(%0.0f)+new(%0.0f)= %0.0f\n',...
                include_file,lr,lc,lr+lc);
            
        end
        
    end

    function do_model_conversions()
        
        % replace ++ -- etc.
        just_slick=@replace_redundancies; %#ok<NASGU>
        scenarios='\+\s*\+|\+\s*\-|\-\s*\+';
        model_eqtns=regexprep(model_eqtns,['\<(\w+\s*)(',scenarios,')(\s*\<\w+)\>'],...
            '${just_slick($1,$2,$3)}');
        
        % replace STEADY_STATE with steady_state
        %----------------------------------------
        model_eqtns=regexprep(model_eqtns,'STEADY_STATE','steady_state');
        
        % replace steady_state(a-b) with (steady_state(a)-b)
        %---------------------------------------------------
        model_eqtns=regexprep(model_eqtns,...
            'steady_state\(\s*(\w+)\s*(\+|\-)\s*(\w+)\s*\)',...
            '(steady_state($1)$2$3)');
        
        if detrend
            % detrend the model
            %------------------
            detrending(deflators)
            
            detrending(log_deflators)
        else
            % replace steady states of log and level variables
            %-------------------------------------------------
            replace_steady_states_for_trending_variables(log_variables,true)
            
            replace_steady_states_for_trending_variables(level_variables,false)
            
        end
        
        % replace ln with log
        model_eqtns=regexprep(model_eqtns,'\<ln\s*\(','log(');
        
        % take care of predetermined variables
        do_predetermined();
        
        % augment model with trend_var and log_trend_var
        %------------------------------
        augment_model_with_trend_variables(trend_var)
        
        augment_model_with_trend_variables(log_trend_var)
        
        % remove or comment out equation tags
        %------------------------------------
        process_equation_tags()
        
        % push standard deviations into equations directly
        %-------------------------------------------------
        add_stdev_to_model()
        
        function detrending(batch)
            
            vnames=fieldnames(batch);
            
            if isempty(vnames)
                
                return
                
            end
            
            c2m_endo_list=parser.cell2matize(endo_list(:).');
            
            switch inputname(1)
                
                case 'deflators'
                    
                    myreplace=@engine_deflators; %#ok<NASGU>
                    
                case 'log_deflators'
                    
                    myreplace=@engine_log_deflators; %#ok<NASGU>
                    
                otherwise
                    
                    error('This should never happen')
                    
            end
            
            cvlist=parser.cell2matize(vnames(:).');
            
            model_eqtns=regexprep(model_eqtns,['(?<!steady_state\s*\(\s*)\<',...
                cvlist,'\>\s*(\([\+\-0-9]+\))?'],'${myreplace($1,$2)}');
            
            function out=engine_deflators(v,lagLead)
                
                [specific_deflator,append]=relead_relag(batch.(v),lagLead);
                
                out=['(',v,append,'*(',specific_deflator,'))'];
                
            end
            
            function out=engine_log_deflators(v,lagLead)
                
                [specific_deflator,append]=relead_relag(batch.(v),lagLead);
                
                out=['(',v,append,'+(',specific_deflator,'))'];
                
            end
            
            function [newdefl,append]=relead_relag(defl,lagLead)
                
                newdefl=defl;
                
                append='';
                
                if isempty(lagLead),return,end
                
                lagLead(isspace(lagLead))=[];
                
                lagLead=strrep(lagLead,'(','');
                
                lagLead=strrep(lagLead,')','');
                
                if strcmp(lagLead,'0'),return,end
                
                append=['{',lagLead,'}'];
                
                newdefl=regexprep(newdefl,['\<',c2m_endo_list,'\>'],['$1',append]);
                
            end
            
        end
        
        function replace_steady_states_for_trending_variables(list,islogvar)
            
            if ~isempty(list)
                % when the model is detrended, the log_variables only contain
                % trend_vars, which do not appear in the model with their
                % steady states. So, no instances of the regexp will be found.
                
                lvv=parser.cell2matize(list(:).');
                
                if islogvar
                    
                    myengine=@myreplace_log_vars; %#ok<NASGU>
                    
                else
                    
                    myengine=@myreplace_level_vars; %#ok<NASGU>
                    
                end
                
                model_eqtns=regexprep(model_eqtns,['steady_state\(\s*',lvv,'\s*\)'],'${myengine($1)}');
                
            end
            
            function out=myreplace_log_vars(v)
                
                out=['steady_state(',v,')'];
                
                out=['(',out,'/(',deflators.(v),'))'];
                
            end
            
            function out=myreplace_level_vars(v)
                
                out=['steady_state(',v,')'];
                
                out=['(',out,'-(',log_deflators.(v),'))'];
                
            end
            
        end
        
        function out=replace_redundancies(a,b,c)
            
            warning(['"',a,b,c,'" not tidy!!!'])
            
            a(isspace(a))=[];
            
            b(isspace(b))=[];
            
            c(isspace(c))=[];
            
            switch b
                
                case '++'
                    
                    b='+';
                    
                case {'-+','+-'}
                    
                    b='-';
                    
                otherwise 
                    
                    error('junior.maih@gmail.com messed up again, please let him know')
            
            end
            
            out=[a,b,c];
            
        end
        
        function do_predetermined()
            
            if ~isempty(pred_vars)
                
                just_neat=@myreplace; %#ok<NASGU>
                
                % add {0} to guys that are current i.e. not followed by a (
                %----------------------------------------------------------
                xpress0=parser.cell2matize(pred_vars);
                xpress=['\<',xpress0,'\>(?!\()'];
                model_eqtns=regexprep(model_eqtns,xpress,'$1{0}');
                
                % Turn all guys with (+-d) into {+-d}
                %------------------------------------
                xpress=['\<',xpress0,'\>\((\+|\-)?(\d+)\)'];
                model_eqtns=regexprep(model_eqtns,xpress,'$1{$2$3}');
                
                % Now substract 1 to the time indices
                %-------------------------------------
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
        
        function process_equation_tags()
            
            % comment out equation tags and push the equation to the next line
            %-----------------------------------------------------------------
            
            tags=regexp(model_eqtns,'[\s*name\s*=','start');
            
            ntags=numel(tags);
            
            for it=ntags:-1:1
                
                tt=tags(it);
                
                nopen=1;
                
                for icol=tt+1:length(model_eqtns)
                    
                    if model_eqtns(icol)=='['
                        
                        nopen=nopen+1;
                        
                    elseif model_eqtns(icol)==']'
                        
                        nopen=nopen-1;
                        
                    end
                    
                    if nopen==0
                        
                        break
                        
                    end
                    
                end
                
                model_eqtns=[model_eqtns(1:tt-1),'%',...
                    model_eqtns(tt:icol),eol, model_eqtns(icol+1:end)];
                
            end
            
        end
        
        function augment_model_with_trend_variables(tv_or_log_tv)
            
            newvars=fieldnames(tv_or_log_tv);
            
            if isempty(newvars)
                
                return
                
            end
            
            if strcmp(inputname(1),'trend_var')
                
                separator='/';
                
            else
                
                separator='-';
                
            end
            
            newvars__=cell2mat(strcat(newvars(:).',', '));
            
            % replace trigger
            trigstart=strfind(endo_block,'endogenous');
            
            trigstart=trigstart(1);
            
            trigend=trigstart+length('endogenous')-1;
                        
            endo_block=[endo_block(trigstart:trigend),' ',newvars__,endo_block(trigend+1:end)];
            
            neweqtns='';
            
            for iii=1:numel(newvars)
                
                v=newvars{iii};
                
                eqtn=[v,separator,v,'{-1}=',tv_or_log_tv.(v),';'];
                
                neweqtns=sprintf('%s\n%s',neweqtns,eqtn);
                
            end
            
            model_eqtns=sprintf('%s\n%s',neweqtns,model_eqtns);
            
        end
        
        function add_stdev_to_model()
            
            expre=cell2mat(strcat(exo_list,'|'));
            
            expre=['\<(',expre(1:end-1),')\>'];
            
            repl=[stderr_name,'_$1*$1'];
            
            model_eqtns=regexprep(model_eqtns,expre,repl);
            
        end
        
    end

end

function raw_code = read_file(dynFileName,my_contains)

if ~my_contains(dynFileName,'.')
    
    success=false;
    
    xtens={'.dyn','.mod'};
    
    for ix=1:numel(xtens)
        
        tmp=[dynFileName,xtens{ix}];
        
        if exist(tmp,'file')
            
            success=true;
            
            dynFileName=tmp;
            
            break
            
        end
        
    end
    
    if ~success
        
        error('No file "%s" with extension .dyn or .mod found',dynFileName)
        
    end
    
elseif ~exist(dynFileName,'file')
    
    error('Unable to find ''%s''.',dynFileName);
    
end

fid = fopen(dynFileName,'r');

if fid == -1
    
    error('Unable to open ''%s'' for reading.',dynFileName);
    
end

raw_code = transpose(fread(fid,'*char'));

fclose(fid);

% remove block comments
%-----------------------
lraw=length(raw_code);

raw_code=stepwise_removal_of_block_comments(raw_code);

lrise=length(raw_code);

fprintf(1,'Removing block comments in "%s": Old(%0.0f), New(%0.0f)\n',dynFileName,lraw,lrise);

% Remove line comments.
%----------------------
lrise_old=lrise;

raw_code = regexprep(raw_code,'(//|%)(.*?\n)','\n');

lrise=length(raw_code);

fprintf(1,'Removing comment lines in "%s": Old(%0.0f), New(%0.0f)\n',dynFileName,lrise_old,lrise);


    function raw_code=stepwise_removal_of_block_comments(raw_code)
        
        % double slash with /* is problematic...
        raw_code=strrep(raw_code,'//','%');
        
        starts=strfind(raw_code,'/*'); %starts=regexp(raw_code,'^/\*','start'); % starts=regexp(raw_code,'(?<!/)/\*','start');
                
        trim_starts()
        
        finishes=strfind(raw_code,'*/');%finishes=regexp(raw_code,'^\*/','start'); %
        
        if numel(starts)~=numel(finishes)
            
            error('mismatches of block comments')
            
        end
        
        while ~ isempty(starts)
            
            blk=starts(end):finishes(end)+1;
            
            if isempty(blk)
                
                raw_code(finishes(end):starts(end))
                
                error('finishing occurring before opening in the code above')
                
            end
            
            raw_code(blk)=[];
            
            starts(end)=[];
            
            finishes(end)=[];
            
        end
        
        function trim_starts()
            
            nstarts=numel(starts);
            
            bad=false(1,nstarts);
            
            for ii=1:nstarts
                
                ping=starts(ii)-1;
                
                if ping>0
                    
                    before=raw_code(ping);
                    
                    bad(ii)=~isspace(before);
                    
                end
                
            end
            
            starts(bad)=[];
            
        end
        
    end


end

function est=extract_estimation(estim_block,std_type,obs_list)

% remove rows starting with corr
% corr_rows=regexp(estim_block,'corr[^;]+;','match');

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

function m=untag_model(m)

taglocs=regexp(m,'%?\[\s*name\s*=','start'); % ?= 0 or 1 time

ntags=numel(taglocs);

for ii=ntags:-1:1
    
    tt=taglocs(ii);
    
    n=1;
    
    % skip "%[" or "["
    if m(tt)=='%'
        
        start=tt+2;
        
    else
        
        start=tt+1;
        
    end
    
    for icol=start:length(m)
        
        if m(icol)=='['
            
            n=n+1;
            
        elseif m(icol)==']'
            
            n=n-1;
            
        end
        
        if n==0
            
            break
            
        end
        
    end
    
    m=[m(1:tt-1),m(icol+1:end)];
    
end

end

function [uList,dupList]=set_list(list)

uList=unique(list);

n=numel(uList);

isdup=false(1,n);

if numel(list)~=n
    
    for ii=1:n
        
        v=uList{ii};
        
        isdup(ii)=sum(strcmp(v,list))>1;
        
    end
    
end

dupList=uList(isdup);

if isempty(uList)
    
    uList={};
    
end

end