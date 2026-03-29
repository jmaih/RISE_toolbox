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
l=0.5*(1-prob);
r=1-l;
for ii=1:npars
    d=Params{ii,2};
    if strcmp(d,'inv_gamma'),d='igamma';end
    o_engine=rdist.(d)(Params{ii,3},Params{ii,4});
    bounds(ii,:)=o_engine([l,r],'icdf');
end
disp(bounds)
%% recover means and standard deviations
clc
MOMS=nan(npars,2);
MOMS2=nan(npars,2);
for ii=1:npars
    d=Params{ii,2};
    if strcmp(d,'inv_gamma'),d='igamma';end
    o_engine=rdist.(d)(bounds(ii,1),bounds(ii,2),prob);
    MOMS(ii,1)=o_engine([],'mean');
    MOMS(ii,2)=o_engine([],'sd');
end
%% display and compare
disp('================ Means and Standard Deviations ================ ')
disp(num2cell(MOMS))
disp('================ Discrepancies ================')
disp(abs(cell2mat(Params(:,3:end))-MOMS))

