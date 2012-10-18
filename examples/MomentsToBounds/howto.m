%% computing the bounds from the mean and standard deviations from dynare
prob=0.9;
Params={
	'alp'  ,  'beta'     ,    0.356000, 0.02	  
	'bet'  ,  'beta'     ,    0.993000, 0.002		
	'gam'  ,  'normal'   ,    0.008500, 0.003	
	'mst'  ,  'normal'   ,    1.000200, 0.007	
	'rho'  ,  'beta'     ,    0.129000, 0.223		
	'psi'  ,  'beta'     ,    0.650000, 0.05		
	'del'  ,  'beta'     ,    0.010000, 0.005		
	'sig_a',  'inv_gamma',    0.035449, inf	
	'sig_m',  'inv_gamma',    0.008862, inf
	};
npars=size(Params,1);
bounds=nan(npars,2);
for ii=1:npars
    bounds(ii,:)=distributions.find_bounds(Params{ii,2},Params{ii,3},Params{ii,4},prob);
end
disp(bounds)
%% recover means and standard deviations
MOMS=nan(npars,2);
for ii=1:npars
    [a,b,moments,fval]=distributions.(Params{ii,2})(bounds(ii,1),bounds(ii,2),prob);
    MOMS(ii,1)=moments.mean;
    MOMS(ii,2)=moments.sd;
end
%% display and compare
disp('================ Means and Standard Deviations ================ ')
disp(num2cell(MOMS))
disp('================ Discrepancies ================')
disp(abs(cell2mat(Params(:,3:end))-MOMS))

