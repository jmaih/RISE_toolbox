function [dictionary,...
    dyn,...
    stat,...
    defs,...
    shadow_tvp,...
    shadow_complementarity]=shadowize(dictionary,AllModels,equation_type)

% dynamic to be merged later with the other portion
%---------------------------------------------------
dyn=struct();
dyn.model=cell(0,1);
dyn.shadow_model=cell(0,1);

% static to be merged later with the other portion
%-------------------------------------------------
stat=struct();
stat.model=cell(0,1);
stat.shadow_model=cell(0,1);
stat.shadow_steady_state_model=cell(0,1);

% plan_syst to be merged later with the other portion in
% dictionary.planner_system 
%---------------------------------------------------------------
shadow_plan_syst=cell(0,1);

% definitions
%-------------
defs=struct();
defs.shadow=cell(0,1);
defs.original=cell(0,1);

shadow_tvp=cell(0,1);
shadow_complementarity=cell(0,1);

for ii=1:numel(equation_type)
    % we don't need the semicolons anymore
    eq_i=AllModels{ii,1};
    o_m='';s_m='';sh_o=''; sh_s='';
    sh_d='';sh_tvp='';sh_ssm='';
    sh_pl='';o_d='';
    sh_mcp='';
    is_def=equation_type(ii)==2;
    is_tvp=equation_type(ii)==3;
    is_mcp=equation_type(ii)==4;
    is_sseq=any(equation_type(ii)==[5,7]);
    is_planner=equation_type(ii)==6;
    
    Vss_pos=0;
    for jj=1:size(eq_i,2)
        item=eq_i{1,jj};
        lead_or_lag=eq_i{2,jj};
        [status,pos]=dictionary.determine_status(item,dictionary);
        switch status
            case 'y'
                if any([is_sseq,is_planner,is_tvp])
                    % is_tvp is added here because only contemporaneous
                    % variables enter the tvp be it in steady state or
                    % during filtering.
                    index=pos;
                elseif is_mcp
                    % this will be used in simulations and therefore
                    % requires the final solution (the final list of
                    % endogenous): the position should be re-assessed
%                     old_pos=pos;
                    pos=find(strcmp(item,{dictionary.endogenous.name}));
                    index=pos;
                else
                    index=abs(lead_or_lag-2);
                    index=dictionary.lead_lag_incidence.before_solve(pos,index);
                end
                if is_def
                    error([mfilename,':: definitions cannot contain variables']);
                elseif is_tvp
                    sh_tvp=[sh_tvp,'y(',sprintf('%0.0f',index),')'];  %#ok<*AGROW>
                    dictionary.endogenous(pos).is_trans_prob=true;
                elseif is_sseq
                    sh_ssm=[sh_ssm,'y(',sprintf('%0.0f',index),')'];
                    sh_ssm=update_time(sh_ssm,lead_or_lag);
                elseif is_planner
                    sh_pl=[sh_pl,'y(',sprintf('%0.0f',index),')'];
                elseif is_mcp
                    sh_mcp=[sh_mcp,'y(',sprintf('%0.0f',index),')'];
                else
                    o_m=[o_m,item];
                    if lead_or_lag
                        o_m=[o_m,'{',sprintf('%0.0f',lead_or_lag),'}'];
                    end
                    sh_o=[sh_o,'y(',sprintf('%0.0f',index),')'];
                    
                    s_m=[s_m,item];
                    sh_s=[sh_s,'y(',sprintf('%0.0f',pos),')'];
                    sh_s=update_time(sh_s,lead_or_lag);
                end
            case 'x'
                if is_def||is_mcp|| is_tvp|| is_planner
                    error([mfilename,':: definitions, complementarity restrictions, endogenous probabilities and planner objective cannot contain shocks. Use auxiliary variables if necessary']);
                elseif is_sseq
                    sh_ssm=[sh_ssm,'x(',sprintf('%0.0f',pos),')'];
                else
                    o_m=[o_m,item];
                    sh_o=[sh_o,'x(',sprintf('%0.0f',pos),')'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,'x(',sprintf('%0.0f',pos),')'];
                end
            case 'param'
                if is_def
                    if lead_or_lag % 'param(',sprintf('%0.0f',pos),')'
                        error('leads on parameters not allowed in definitions')
                    end
                    sh_d=[sh_d,sprintf('param(%0.0f)',pos)];%
                    o_d=[o_d,item];
                elseif is_tvp
                    if lead_or_lag
                        error('leads on parameters not allowed in time-varying probabilities')
                    end
                    sh_tvp=[sh_tvp,sprintf('param(%0.0f)',pos)];
                elseif is_sseq
                    if lead_or_lag
                        error('leads on parameters not allowed in steady state equations')
                    end
                    sh_ssm=[sh_ssm,sprintf('param(%0.0f)',pos)];
                elseif is_planner
                    if lead_or_lag
                        error('leads on parameters not allowed in the planner objective')
                    end
                    sh_pl=[sh_pl,sprintf('param(%0.0f)',pos)];
                elseif is_mcp % leads and lags already checked at capture
                    sh_mcp=[sh_mcp,sprintf('param(%0.0f)',pos)];
                else
                    pname='param';
                    if lead_or_lag
                        pname='sparam';
                        if ~dictionary.parameters(pos).is_switching
                            error(['"',dictionary.parameters(pos).name,'" declared as a switching parameter but found to have a lead'])
                        end
                    end
                    o_m=[o_m,item];
                    sh_o=[sh_o,pname,'(',sprintf('%0.0f',pos),')'];
                    s_m=[s_m,item];
                    sh_s=[sh_s,sprintf('param(%0.0f)',pos)];
                end
            case 'def'
                if is_def
                    sh_d=[sh_d,'def(',sprintf('%0.0f',pos),')'];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,'def(',sprintf('%0.0f',pos),')'];
                elseif is_sseq
                    sh_ssm=[sh_ssm,'def(',sprintf('%0.0f',pos),')'];
                elseif is_planner
                    sh_pl=[sh_pl,'def(',sprintf('%0.0f',pos),')'];
                elseif is_mcp
                    sh_mcp=[sh_mcp,'def(',sprintf('%0.0f',pos),')'];
                else
                    o_m=[o_m,item];
                    s_m=[s_m,item];
                    sh_o=[sh_o,'def(',sprintf('%0.0f',pos),')'];
                    sh_s=[sh_s,'def(',sprintf('%0.0f',pos),')'];
                end
            case 'tvp'
                if ~is_tvp
                    error([mfilename,':: endogenous probability cannot be used elsewhere than in their block'])
                end
                sh_tvp=[sh_tvp,item];
            otherwise
                if strcmp(item,'steady_state')
                    % then two blocks later is the name of the variable
                    Vss_pos=jj+2;
                    Vss=eq_i{1,Vss_pos};
                    if is_mcp
                        pos=find(strcmp(Vss,{dictionary.endogenous.name}));
                    else
                        pos=find(strcmp(Vss,{dictionary.endogenous.name}));
                    end
                    eq_i{1,jj+2}=sprintf('%0.0f',pos);
                end
                if is_def
                    sh_d=[sh_d,item];
                    o_d=[o_d,item];
                elseif is_tvp
                    sh_tvp=[sh_tvp,item];
                elseif is_sseq
                    sh_ssm=[sh_ssm,item];
                elseif is_planner
                    sh_pl=[sh_pl,item];
                elseif is_mcp
                    sh_mcp=[sh_mcp,item];
                else
                    if isequal(jj,Vss_pos)
                        o_m=[o_m,Vss];
                        s_m=[s_m,Vss];
                    else
                        o_m=[o_m,item];
                        s_m=[s_m,item];
                    end
                    sh_o=[sh_o,item];
                    sh_s=[sh_s,item];
                end
        end
    end
    if is_def
        defs.shadow=[defs.shadow;{sh_d}]; % put back the semicolon as this is going to be evaluated
        defs.original=[defs.original;{o_d}];
    elseif is_tvp
        shadow_tvp=[shadow_tvp;{sh_tvp}];
    elseif is_sseq
        stat.shadow_steady_state_model=[stat.shadow_steady_state_model;{sh_ssm}];
        % put back the semicolon as this is going to be evaluated
    elseif is_planner
        shadow_plan_syst=[shadow_plan_syst;{sh_pl}];
    elseif is_mcp
        shadow_complementarity=[shadow_complementarity;{sh_mcp}];
    else
        dyn.model=[dyn.model;{o_m}];
        dyn.shadow_model=[dyn.shadow_model;{sh_o}];
        stat.model=[stat.model;{s_m}];
        stat.shadow_model=[stat.shadow_model;{sh_s}];
    end
end

if isfield(dictionary.planner_system,'shadow_model')
    if ~isempty(shadow_plan_syst)
        error('several planner objectives?')
    end
else
    dictionary.planner_system.shadow_model=shadow_plan_syst;
end

% replace the steady state calls
%-------------------------------
dyn.shadow_model=parser.replace_steady_state_call(dyn.shadow_model);
stat.shadow_model=parser.replace_steady_state_call(stat.shadow_model);
shadow_complementarity=parser.replace_steady_state_call(shadow_complementarity);
% remove the useless call to the steady state in steady state models
%---------------------------------------------------------------------
stat.shadow_model=strrep(stat.shadow_model,'ss','y');

end

function sh=update_time(sh,lead_or_lag)
switch lead_or_lag
    case -1
        sh=[sh,'{-1}'];
    case 1
        sh=[sh,'{+1}'];
end
end

