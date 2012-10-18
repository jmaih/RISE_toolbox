% gives approximate values for the steady state of fs2000 
% this is the equivalent of dynare's initvals
function [ys,retcode]=fs2000_steadystate_initval(~,flag)

retcode=0;

% var_names={'P','R','W','c','d','dA','e','gp_obs','gy_obs','k','l','m','n','y'};
% 
% ys =[2.2582,1.0212,4.5959,0.4477,0.8494,1.0030,1,1.0080,1.0030,5.8012,0.8604,1.0110,0.1872,0.5808]';
%     
tmp={
	'P'			,	    2.258154387910923
%     'P_AUX_F1'	,	    2.258154387910923
    'R'			,	    1.021212121212121
    'W'			,	    4.595903784741778
    'c'			,	    0.447710752379204
%     'c_AUX_F1'	,	    0.447710752379204
    'd'			,	    0.849424911502341
    'dA'		,	    1.003004504503377
    'e'			,	    1.000000000000000
    'gp_obs'	,	    1.007971544953910
    'gy_obs'	,	    1.003004504503377
    'k'			,	    5.801216036088844
    'l'			,	    0.860424911502341
    'm'			,	    1.011000000000000
    'n'			,	    0.187215605852959
    'y'			,	    0.580765090448550
};

switch flag
    case 0
        ys=tmp(:,1);
    case 1
        ys =cell2mat(tmp(:,2));
    otherwise
end

    
