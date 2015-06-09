function varargout=load_data(obj)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

just_starting=isempty(obj);
[varargout{1:nargout}]=load_data@rise_generic(obj);

if just_starting
    return
end

obj=varargout{1};
data_provided=obj.options.data.NumberOfVariables>0;
if data_provided %simulation_available || 
    estim_start_date=obj.options.estim_start_date;
    obj.dates_filtering=serial2date(date2serial(estim_start_date)+(0:obj.data.nobs));
    obj.dates_smoothing=obj.dates_filtering(1:end-1);

    if obj.is_dsge_var_model
        data=obj.data.y(:,obj.data.start,obj.data.finish);
        const=obj.options.dsge_var_constant;
        n=obj.NumberOfObservables(1); % endogenous observables
        p=obj.options.dsge_varlag;
        Y=data(p+1:end,:);
        smpl=size(Y,1);
        
        X=nan(smpl,const+n*p);
        if const
            X(:,1)=1;
        end
        for ii=1:p
            X(:,const+((ii-1)*n+1:ii*n))=data(p+1-ii:end-ii,:);
        end
        dsge_var=struct('YY',Y'*Y,'YX',Y'*X,...
            'XX',X'*X,'XY',X'*Y,'T',smpl,...
            'n',n,'p',p,'constant',const);
        obj=obj.set_properties('dsge_var',dsge_var);
    end
     
    % information on the conditional variables
    %------------------------------------------
    kmax=max(obj.exogenous.shock_horizon);
    if kmax
        if obj.data.npages<kmax+2
            warning('the anticipation horizon of agents will be reduced since it exceeds the number of advance information')
        end
        obj.exogenous.shock_horizon=min(obj.exogenous.shock_horizon,kmax+2);
    end
end

varargout{1}=obj;

end