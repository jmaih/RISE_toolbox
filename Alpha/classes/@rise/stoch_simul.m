function [oo_]=stoch_simul(obj,var_list,varargin)%,omega
if isempty(obj)
    oo_=struct('stoch_sim_hpfilter_lambda',false);
    return
end
oo_=struct();

if nargin<2
    var_list=[];
end
% should also allow for passing options
% get that small set of original variables

test_for_deep_parameters_calibration();

obj=set_options(obj,varargin{:});
if obj.is_linear_model
    obj=set_options(obj,'solve_order',1);
end
if obj.options.solve_order == 1
    obj=set_options(obj,'irf_draws',1);% replic
end

if isempty(var_list)
    var_list=get(obj,'endo_list(original)');
end

% iter_ = max(obj.options.periods,1);
% if obj.exogenous.number > 0
%     oo_.exo_simul= ones(iter_ + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
% end

Lead_lag_incidence=obj.Lead_lag_incidence;
check_model();

[obj,retcode] = solve(obj);

if retcode
    decipher_error(retcode);
    return
end

% if ~obj.options.noprint
% Model summary
%--------------
disp(' ')
disp('MODEL SUMMARY')
disp(' ')
disp(['  Number of variables:         ' int2str(obj.endogenous.number(2))])
disp(['  Number of stochastic shocks: ' int2str(sum(obj.exogenous.number))])

disp(['  Number of state variables:   ' int2str(nnz(Lead_lag_incidence(:,3)))])
disp(['  Number of jumpers:           ' int2str(nnz(Lead_lag_incidence(:,1)))])

static=sum(Lead_lag_incidence(:,3)==0 & Lead_lag_incidence(:,1)==0);
disp(['  Number of static variables:  ' int2str(static)])

disp(' ')
if obj.options.solve_order <= 2
    print_solution(obj,var_list,true);
end
% end


% if obj.options.irf
oo_.irfs=irf(obj);
% end

if obj.options.simul_periods > 0
    db= simulate(obj);
    disp_moments(db,var_list);
    oo_.simulations=db;
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

    function check_model()
        xlen = 0;% DynareModel.maximum_exo_lag+DynareModel.maximum_exo_lead + 1;
        if ~ all(Lead_lag_incidence(:,2)) > 0
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

    function disp_moments(db,var_list)
        if nargin<2
            var_list=[];
        end
        if isempty(var_list)
            var_list=get(obj,'endo_list(original)');
        end
        stoch_sim_hpfilter_lambda=obj.options.stoch_sim_hpfilter_lambda;
        
        ivar=locate_variables(var_list,db.varnames);
        oo_.mean=mean(db);
        
        % warning_old_state = warning;
        % warning off
        %
        
        if stoch_sim_hpfilter_lambda
            [oo_.hp_trend,oo_.hp_cycle] = hpfilter(db,stoch_sim_hpfilter_lambda);
            db=oo_.hp_cycle;
        else
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
        if stoch_sim_hpfilter_lambda
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
B=concatenate(data,'%8.6f');
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
