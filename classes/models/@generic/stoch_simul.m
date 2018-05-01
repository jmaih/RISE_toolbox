function [oo_]=stoch_simul(obj,var_list,varargin)%,omega
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

if isempty(obj)
    
    oo_=cell(0,4);
    
    return
    
end

if nargin<2
    
    var_list=[];
    
end

nobj=numel(obj);

if nobj>1
    
    oo_=cell(1,nobj);
    
    for iobj=1:nobj
        
        fprintf('\n ******************* model(%0.0f): %s ******************* \n',iobj,obj(iobj).filename)
        
        [oo_{iobj}]=stoch_simul(obj(iobj),var_list,varargin{:});
        
    end
    
    return
    
end

oo_=struct();
% should also allow for passing options
% get that small set of original variables

test_for_deep_parameters_calibration();

obj=set(obj,varargin{:});

if isa(obj,'dsge')
    
    if obj.options.solve_order == 1
        
        obj=set(obj,'irf_draws',1);% replic
        
    end
    
end

if isempty(var_list)
    
    if isa(obj,'dsge')
        
        var_list=get(obj,'endo_list(original)');
        
        lead_lag_incidence=obj.lead_lag_incidence;
        
        check_model();
        
    else
        
        var_list=get(obj,'endo_list');
        
    end
    
end

% iter_ = max(obj.options.periods,1);
% if obj.exogenous.number > 0
%     oo_.exo_simul= ones(iter_ + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
% end

[obj,retcode] = solve(obj);

if retcode
    
    utils.error.decipher(retcode);
    
    return
    
end

% if ~obj.options.noprint
% Model summary
%--------------
disp(' ')
disp('MODEL SUMMARY')
disp(' ')
disp(['  Number of variables:         ' int2str(obj.endogenous.number)])
disp(['  Number of stochastic shocks: ' int2str(obj.exogenous.number(1))])
if isa(obj,'dsge')
disp(['  Number of state variables:   ' int2str(nnz(lead_lag_incidence.after_solve(:,3)))])
disp(['  Number of jumpers:           ' int2str(nnz(lead_lag_incidence.after_solve(:,1)))])

static=sum(lead_lag_incidence.after_solve(:,3)==0 & lead_lag_incidence.after_solve(:,1)==0);
disp(['  Number of static variables:  ' int2str(static)])
end
disp(' ')
if isa(obj,'dsge') && obj.options.solve_order <= 2
    print_solution(obj,var_list,[],true);
end
% end


% if obj.options.irf
oo_.irfs=irf(obj);
% end

if obj.options.simul_periods > 0
    db= simulate(obj);
    db_names=disp_moments(db);
    oo_.simulations=db;
    do_storage()
else
    for ii=1:10
        disp('No theoretical moments implemented yet. Please send a reminder to junior.maih@gmail.com')
    end
    return
    %     disp_th_moments(var_list);
end

% if obj.options.SpectralDensity.trigger == 1
%     [omega,f] = UnivariateSpectralDensity(oo_.dr,var_list);
% end

    function do_storage()
        fields=fieldnames(oo_);
        endo_names=obj.endogenous.name;
        pos=locate_variables(endo_names,db_names);
        for ifield=1:numel(fields)
            name=fields{ifield};
            if isstruct(oo_.(name))
                continue
            elseif iscell(oo_.(name))
                % this has to be checked first
                for icel=1:numel(oo_.(name))
                    oo_.(name){icel}=oo_.(name){icel}(pos,pos);
                end
            elseif isvector(oo_.(name))
                tmp=struct();
                for iname=1:numel(endo_names)
                    tmp.(endo_names{iname})=oo_.(name)(pos(iname));
                end
                oo_.(name)=tmp;
            elseif ismatrix(oo_.(name))
                oo_.(name)=oo_.(name)(pos,pos);
            end
        end
    end

    function check_model()
        xlen = 0;% DynareModel.maximum_exo_lag+DynareModel.maximum_exo_lead + 1;
        if ~ all(lead_lag_incidence.after_solve(:,2)) > 0
            error ('Error in model specification: some variables do not appear as current') ;
        end
        
        if xlen > 1
            error (['stochastic exogenous variables must appear only at the' ...
                ' current period. Use additional endogenous variables']) ;
        end
    end

    function test_for_deep_parameters_calibration()
        nan_params=isnan(obj);
        if ~isempty(nan_params)
            tmp = dbstack;
            disp(' ')
            disp(nan_params)
            message = [tmp(2).name,':: The parameters above have no value. ',...
                ' If they are not initialized in a steadystate file, RISE ',...
                ' may not be able to solve the model' ];
            message_id  = 'RISE:ParameterCalibration:NaNValues';
            warning(message_id,message);
        end
    end

    function db_names=disp_moments(db)
        
        if isempty(var_list)
            
            var_list=get(obj,'endo_list(original)');
            
        end
        
        stoch_sim_hpfilter_lambda=obj.options.simul_hpfilter_lambda;
        
        if ~isa(db,'ts')
            
            db=ts.collect(db);
            
        end
        
        db_names=db.varnames;
        ivar=locate_variables(var_list,db_names);
        oo_.mean=mean(db);
        
        % warning_old_state = warning;
        % warning off
        %
        
        if isempty(stoch_sim_hpfilter_lambda)
            
            db = bsxfun(db,@minus,oo_.mean);
            
        end
        oo_.vcov=cov(db);
        oo_.skewness=skewness(db);
        oo_.kurtosis=kurtosis(db);
        oo_.variance = diag(oo_.vcov);oo_.variance=oo_.variance(:)';
        oo_.stdev = sqrt(oo_.variance);
        oo_.corrcoef = corrcoef(db);
        
        % if options_.nomoments == 0
        title='MOMENTS OF SIMULATED VARIABLES';
        
        if ~isempty(stoch_sim_hpfilter_lambda)
            
            title = [title,'(HP filter, lambda = ',num2str(stoch_sim_hpfilter_lambda),')'];
            
        end
        
        data=[{'VARIABLE','MEAN','STD. DEV.','VARIANCE','SKEWNESS','KURTOSIS'}
            [var_list(:),num2cell([ oo_.mean(ivar)',oo_.stdev(ivar)',...
            oo_.variance(ivar)',oo_.skewness(ivar)',oo_.kurtosis(ivar)'])]
            ];
        
        reprint(data,title);
        
        % if options_.nocorr == 0
        %     if options_.noprint == 0
        title = 'CORRELATION OF SIMULATED VARIABLES';
        
        if stoch_sim_hpfilter_lambda
            title = [title ' (HP filter, lambda = ' ...
                num2str(stoch_sim_hpfilter_lambda) ')'];
            
        end
        
        data=[[{'VARIABLE'},var_list(:)']
            [var_list(:),num2cell(oo_.corrcoef(ivar,ivar))]
            ];
        reprint(data,title);
        
        %     end
        % end
        
        ar = obj.options.autocorr_ar;
        if ar > 0
            autocorr = [];
            y=double(db);
            for i=1:ar
                oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
                autocorr = [ autocorr diag(oo_.autocorr{i}) ];
            end
            %     if options_.noprint == 0
            title = 'AUTOCORRELATION OF SIMULATED VARIABLES';
            if stoch_sim_hpfilter_lambda
                title = [title,'(HP filter, lambda = ',num2str(stoch_sim_hpfilter_lambda),')'];
            end
            data=[[{'Order'},num2cell(1:ar)]
                [var_list(:),num2cell(autocorr(ivar,:))]
                ];
            reprint(data,title);
            %     end
        end
        
        % warning(warning_old_state);
    end
end
function reprint(data,title)
B=utils.table.concatenate(data,'%8.6f');
body_format='';
% start from the end
for bb=size(B,2):-1:1
    body_format=['%',int2str(B{2,bb}),'s ',body_format]; %#ok<AGROW>
end
nrows=size(B{1,1},1);
number_of_headers=size(B,2);
mycell={sprintf('\n%s ',title)};
for rr=1:nrows
    data_ii=cell(1,number_of_headers);
    for jj=1:number_of_headers
        data_ii{jj}=B{1,jj}(rr,:);
    end
    mycell=[mycell;{sprintf(body_format,data_ii{:})}];
end
for irow=1:numel(mycell)
    fprintf(1,'%s \n',mycell{irow});
end
end
