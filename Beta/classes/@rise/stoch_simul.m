function [obj,irfs,omega]=stoch_simul(obj,var_list)

if isempty(obj)
    obj=struct();
    return
end

% should also allow for passing options
% get that small set of original variables
 
Lead_lag_incidence=obj.Lead_lag_incidence;
disp(' ')
disp('MODEL SUMMARY')
disp(' ')
disp(['  Number of variables:         ' int2str(obj.NumberOfEndogenous(2))])
disp(['  Number of stochastic shocks: ' int2str(obj.NumberOfExogenous)])

disp(['  Number of state variables:   ' int2str(nnz(Lead_lag_incidence(:,3)))])
disp(['  Number of jumpers:           ' int2str(nnz(Lead_lag_incidence(:,1)))])

static=Lead_lag_incidence(:,3)==0 & Lead_lag_incidence(:,1)==0;
disp(['  Number of static variables:  ' int2str(static)])

disp('MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS')
varexo={obj.varexo.name};
disp('FIX THE ISSUE WITH THE OBSERVED EXOGENOUS HERE: GIVE THEM ZERO VARIANCE')
disp([{},varexo;varexo',num2cell(eye(obj.NumberOfExogenous))])
disp(' ')
obj.print_solution(var_list);

obj=obj.simulate;
 
if options_.nomoments == 0
   disp_moments(oo_.endo_simul,var_list);
end


if options_.irf

end

if options_.SpectralDensity.trigger == 1
    [omega,f] = UnivariateSpectralDensity(oo_.dr,var_list);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = [ m' s' s2' (mean(y.^3)./s2.^1.5)' (mean(y.^4)./(s2.*s2)-3)' ];    
title='MOMENTS OF SIMULATED VARIABLES';
headers={'VARIABLE','MEAN','STD. DEV.','VARIANCE','SKEWNESS','KURTOSIS'};
disp([headers;var_list',num2cell(z)])

corr = (y'*y/size(y,1))./(s'*s);
disp('CORRELATION OF SIMULATED VARIABLES');
headers = {'VARIABLE',var_list};
disp([headers;var_list',num2cell(corr)])

	 
 ar = options_.ar;
 if ar > 0
     autocorr = [];
     for i=1:ar
         oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
         autocorr = [ autocorr diag(oo_.autocorr{i}) ];
     end
     if options_.noprint == 0
         title = 'AUTOCORRELATION OF SIMULATED VARIABLES';
         headers = char('VARIABLE',int2str([1:ar]'));
         dyntable(title,headers,labels,autocorr,size(labels,2)+2,8,4);
     end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [A,info]=theoretical_autocorrelations(obj,ar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[oo_.gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_,options_);
 m = dr.ys(ivar);
 non_stationary_vars = setdiff(1:length(ivar),stationary_vars);
 m(non_stationary_vars) = NaN;
 
 i1 = find(abs(diag(oo_.gamma_y{1})) > 1e-12);
 s2 = diag(oo_.gamma_y{1});
 sd = sqrt(s2);
 if options_.order == 2
     m = m+oo_.gamma_y{options_.ar+3};
 end
 
 z = [ m sd s2 ];
 oo_.mean = m;
 oo_.var = oo_.gamma_y{1};
 
 if ~options_.noprint %options_.nomoments == 0
     title='THEORETICAL MOMENTS';
     if options_.hp_filter
         title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];
     end
     headers=char('VARIABLE','MEAN','STD. DEV.','VARIANCE');
     labels = deblank(M_.endo_names(ivar,:));
     lh = size(labels,2)+2;
     dyntable(title,headers,labels,z,lh,11,4);
     if M_.exo_nbr > 1 && size(stationary_vars, 1) > 0
         disp(' ')
         title='VARIANCE DECOMPOSITION (in percent)';
         if options_.hp_filter
             title = [title ' (HP filter, lambda = ' ...
                      num2str(options_.hp_filter) ')'];
         end
         headers = M_.exo_names;
         headers(M_.exo_names_orig_ord,:) = headers;
         headers = char(' ',headers);
         lh = size(deblank(M_.endo_names(ivar(stationary_vars),:)),2)+2;
         dyntable(title,headers,deblank(M_.endo_names(ivar(stationary_vars), ...
                                                      :)),100*oo_.gamma_y{options_.ar+2}(stationary_vars,:),lh,8,2);
     end
     
     conditional_variance_steps = options_.conditional_variance_decomposition;
     if length(conditional_variance_steps)
         oo_ = display_conditional_variance_decomposition(conditional_variance_steps,...
                                                          ivar,dr,M_, ...
                                                          options_,oo_);
     end
 end
 
 if length(i1) == 0
     disp(' ')
     disp('All endogenous are constant or non stationary, not displaying correlations and auto-correlations')
     disp(' ')
     return;
 end
 
 if options_.nocorr == 0 && size(stationary_vars, 1) > 0
     corr = oo_.gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
     if ~options_.noprint,
         disp(' ')
         title='MATRIX OF CORRELATIONS';
         if options_.hp_filter
             title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];
         end
         labels = deblank(M_.endo_names(ivar(i1),:));
         headers = char('Variables',labels);
         lh = size(labels,2)+2;
         dyntable(title,headers,labels,corr,lh,8,4);
     end
 end
 if options_.ar > 0 && size(stationary_vars, 1) > 0
     z=[];
     for i=1:options_.ar
         oo_.autocorr{i} = oo_.gamma_y{i+1};
         z(:,i) = diag(oo_.gamma_y{i+1}(i1,i1));
     end
     if ~options_.noprint,      
         disp(' ')    
         title='COEFFICIENTS OF AUTOCORRELATION';      
         if options_.hp_filter        
             title = [title ' (HP filter, lambda = ' num2str(options_.hp_filter) ')'];      
         end      
         labels = deblank(M_.endo_names(ivar(i1),:));      
         headers = char('Order ',int2str([1:options_.ar]'));
         lh = size(labels,2)+2;
         dyntable(title,headers,labels,z,lh,8,4);
     end  
 end
