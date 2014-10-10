function obj=redo_declarations(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


old_exo_names=obj.exogenous.name;
observed_exo_names=get(obj,'exo_list(observed)');
old_exo_names=setdiff(old_exo_names,observed_exo_names);
exo_names=old_exo_names;
old_endo_names=obj.endogenous.name;
endo_names=old_endo_names;
orig_endo_nbr=obj.endogenous.number(end);

param_template=obj.param_template;
sig_loc=find(strcmp(param_template(1,:),'sig'));
omg_loc=find(strcmp(param_template(1,:),'omg'));
param_template_ac=param_template;
param_template_ac(:,[sig_loc,omg_loc])=[];
newparam_templ=param_template_ac;
endo_locs=struct();
endo_locs.y={old_endo_names,[]};
exo_locs=struct();
exo_locs.E={old_exo_names,[]};
param_locs=struct();
if obj.time_varying_parameters(1) % time-varying autoregressive parameters
    % create variables a and shocks ETA_a
    %------------------------------------
    [endo_names,new_names]=expand_variable_list(endo_names,param_template_ac,[],'_');
    endo_locs.a_={new_names,[]};
    [exo_names,new_names]=expand_variable_list(exo_names,param_template_ac,'ETA_');
    exo_locs.ETA_a={new_names,[]};
    
    % create the rho_a and theta_a parameters
    %-----------------------------------------
    [param_template_theta_ac,new_names]=recreate_template(param_template_ac,'theta_');
    param_locs.theta_a={new_names,[]};
    [param_template_rho_ac,new_names]=recreate_template(param_template_ac,'rho_');
    if obj.random_walk_parameters(1)
        % set parameters a=0, rhoa=1
        for icol=1:size(param_template_ac,2)
            bad=isnan(param_template_ac{2,icol});
            param_template_ac{2,icol}(bad)=0;
            bad=isnan(param_template_rho_ac{2,icol});
            param_template_rho_ac{2,icol}(bad)=1;
        end
    else
        param_locs.rho_a={new_names,[]};
    end
    
    newparam_templ=[newparam_templ,param_template_rho_ac,param_template_theta_ac];
end
[~,new_names]=recreate_template(param_template_ac);
param_locs.a={new_names,[]};

param_template_sig=param_template(:,sig_loc);
newparam_templ=[newparam_templ,param_template_sig];
if obj.time_varying_parameters(2)  % time-varying standard deviations
    % create variables sig and shocks UPSIL_sig
    %--------------------------------------------
    [endo_names,new_names]=expand_variable_list(endo_names,param_template_sig,[],'_');
    endo_locs.sig_={new_names,[]};
    [exo_names,new_names]=expand_variable_list(exo_names,param_template_sig,'UPSIL_');
    exo_locs.UPSIL_sig={new_names,[]};
    
    % create the rho_sig and theta_sig parameters
    %--------------------------------------------
    [param_template_theta_sig,new_names]=recreate_template(param_template_sig,'theta_');
    param_locs.theta_sig={new_names,[]};
    [param_template_rho_sig,new_names]=recreate_template(param_template_sig,'rho_');
    if obj.random_walk_parameters(2)
        % set rho_sig=1
        for icol=1:size(param_template_rho_sig,2)
            bad=isnan(param_template_rho_sig{2,icol});
            param_template_rho_sig{2,icol}(bad)=1;
        end
    else
        param_locs.rho_sig={new_names,[]};
    end
    newparam_templ=[newparam_templ,param_template_rho_sig,param_template_theta_sig];
end
[~,new_names]=recreate_template(param_template_sig);
param_locs.sig={new_names,[]};

param_template_omg=param_template(:,omg_loc);
newparam_templ=[newparam_templ,param_template_omg];
if obj.time_varying_parameters(3)  % time-varying correlations
    % create variables omg and shocks CHI_omg
    %--------------------------------------------
    [endo_names,new_names]=expand_variable_list(endo_names,param_template_omg,[],'_');
    endo_locs.omg_={new_names,[]};
    [exo_names,new_names]=expand_variable_list(exo_names,param_template_omg,'CHI_');
    exo_locs.CHI_omg={new_names,[]};
    
    % create the rho_omg and theta_omg parameters
    %--------------------------------------------
    [param_template_theta_omg,new_names]=recreate_template(param_template_omg,'theta_');
    param_locs.theta_omg={new_names,[]};
    [param_template_rho_omg,new_names]=recreate_template(param_template_omg,'rho_');
    if obj.random_walk_parameters(3)
        % set omg=0 and rho_omg=1;
        for icol=1:size(param_template_omg,2)
            bad=isnan(param_template_omg{2,icol});
            param_template_omg{2,icol}(bad)=0;
            bad=isnan(param_template_rho_omg{2,icol});
            param_template_rho_omg{2,icol}(bad)=1;
        end
    else
        param_locs.rho_omg={new_names,[]};
    end
    newparam_templ=[newparam_templ,param_template_rho_omg,param_template_theta_omg];
end
[~,new_names]=recreate_template(param_template_omg);
param_locs.omg={new_names,[]};

% set the final parameter template
%---------------------------------
obj.param_template=newparam_templ;

% create the parameters
%----------------------
new_params=vartools.create_baseline_parameters(obj.param_template);

endo_names=add_description(endo_names,obj.endogenous.name,obj.endogenous.tex_name);
exo_names=add_description(exo_names,obj.exogenous.name,obj.exogenous.tex_name);
obs_names=obj.observables.name;

% markov chains to be modified if we find a way of the estimating this
% model consistently
%---------------------------------------------------------------------
markchains = struct('name','const','number_of_states',1,...
    'is_endogenous',false,'param_list',{new_params},...
    'param_list_tex',{new_params},'is_switching',false,'duration',Inf);
% collect the block exogenous guys
%---------------------------------
block_exogenous=get(obj,'endo_list(block_exogenous)');
% now redo the endogenous, exogenous, observables, and markov chains
%-------------------------------------------------------------------
obj=rise_generic.reset(obj,endo_names,exo_names,obs_names,markchains);
obj.endogenous.is_block_exogenous=ismember(obj.endogenous.name,block_exogenous);
obj.endogenous.is_original=ismember(obj.endogenous.name,old_endo_names);

locs=struct('endogenous',endo_locs,'exogenous',exo_locs,'parameters',param_locs);
ff=fieldnames(locs);
for ii=1:numel(ff)
    ff_=fieldnames(locs.(ff{ii}));
    for jj=1:numel(ff_)
        variables=locs.(ff{ii}).(ff_{jj}){1}(:)';
        locs.(ff{ii}).(ff_{jj}){1}=variables;
        pos=locate_variables(variables,obj.(ff{ii}).name);
        if any(strcmp(ff_{jj},{'a','a_','omg','omg_'}))
            pos=rematch_positions(pos,variables);
        end
        locs.(ff{ii}).(ff_{jj}){2}=pos(:)';
        obj.([ff{ii},'_positions'])=locs.(ff{ii});
    end
end

    function pairs=rematch_positions(pairs,b_names)
        pairs=num2cell(pairs);
        for local_position=1:numel(b_names)
            name=b_names{local_position};
            state_position=pairs{local_position};
            if ~strcmp(name(end),'_')
                % add an underscore for parameters
                name=[name,'_']; %#ok<AGROW>
            end
            underscores=find(name=='_');
            row=str2double(name(underscores(1)+1:underscores(2)-1));
            col=str2double(name(underscores(2)+1:underscores(3)-1));
            if any(strncmp(name,{'a','o'},1))
                mat=1;
                if name(1)=='a'
                    mat=str2double(name(2:underscores(1)-1));
                end
                matrix_position=(mat-1)*orig_endo_nbr^2+(col-1)*orig_endo_nbr+row;
            elseif name(1)=='c'
                matrix_position=obj.nlags*orig_endo_nbr^2+(col-1)*orig_endo_nbr+row;
            else
                error(['"',name,'" appears to be a wrong variable name'])
            end
            % state_pos,local_pos,matrix_pos
            pairs{local_position}=[state_position,local_position,matrix_position];
        end
    end
end

function [names,new_names]=expand_variable_list(names,tmplate,prefix,suffix)
if nargin<4
    suffix=[];
    if nargin<3
        prefix=[];
    end
end
if isempty(prefix)
    prefix='';
end
if isempty(suffix)
    suffix='';
end
siglist=vartools.create_baseline_parameters(tmplate);
new_names=strcat(prefix,siglist,suffix);
names=[names,new_names(:)'];
end

function [tmplate,plist]=recreate_template(tmplate,prefix,suffix)
if nargin<3
    suffix=[];
    if nargin<2
        prefix=[];
    end
end
if isempty(prefix)
    prefix='';
end
if isempty(suffix)
    suffix='';
end

tmplate(1,:)=strcat(prefix,tmplate(1,:),suffix);
if nargout>1
    plist=vartools.create_baseline_parameters(tmplate);
end
end

function xx=add_description(xx,yy,ydesc)
xx=xx(:);
locs=locate_variables(xx,yy,true);
isgood=~isnan(locs);
xx=[xx,xx];
xx(isgood,2)=ydesc(locs(isgood));
xx(:,2)=strcat('"',xx(:,2),'"');
xx(:,2)=strrep(xx(:,2),'""','"');
xx=xx';
xx=xx(:)';
end