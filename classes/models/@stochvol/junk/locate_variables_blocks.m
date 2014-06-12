function [lambda_id,nu_id,obs_id,B_id,eta_id,e_id,pairs]=locate_variables_blocks(obj)
orig_endo_nbr=numel(get(obj,'endo_list(original)'));
if obj.stochastic_volatility
    lambda_names=regexp(obj.endogenous.name,'lambda_\d+','match');
    lambda_names=[lambda_names{:}];
    lambda_id=locate_variables(lambda_names,obj.endogenous.name);
    nu_id=find(strncmp(obj.exogenous.name,'NU_EXO_',7));
else
    lambda_id=[];
    nu_id=[];
end
obs_id=obj.observables.state_id(obj.observables.is_endogenous);
if obj.time_varying_parameters
    b_names=regexp(obj.endogenous.name,'b(a\d+|c)_\d+_\d+_','match');
    b_names=[b_names{:}];
    B_id=locate_variables(b_names,obj.endogenous.name);
    % use the names to reconstruct the B matrix
    %-------------------------------------------
    pairs=nan(size(B_id));
    for iname=1:numel(b_names)
        name=b_names{iname};
        underscores=find(name=='_');
        row=str2double(name(underscores(1)+1:underscores(2)-1));
        col=str2double(name(underscores(2)+1:underscores(3)-1));
        if name(2)=='a'
            mat=str2double(name(3:underscores(1)-1));
            left=(mat-1)*orig_endo_nbr^2+(col-1)*orig_endo_nbr+row;
        elseif name(2)=='c'
            left=obj.nlags*orig_endo_nbr^2+(col-1)*orig_endo_nbr+row;
        else
            error(['"',name,'" appears to be a wrong variable name'])
        end
        pairs(iname)=iname+left*1i;
    end
    pairs=[imag(pairs),real(pairs)];
    % ETA shocks to parameters
    %-------------------------
    eta_names=regexp(obj.exogenous.name,'ETA_b(a\d+|c)_\d+_\d+','match');
    eta_names=[eta_names{:}];
    eta_id=locate_variables(eta_names,obj.exogenous.name);
else
    B_id=[];
    pairs=[];
    eta_id=[];
end
% structural shocks
%------------------
e_id=find(strncmp(obj.exogenous.name,'EXO_',4));
end
