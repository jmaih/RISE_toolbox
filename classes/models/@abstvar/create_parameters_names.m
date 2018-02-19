function p=create_parameters_names(nvars,lags,constant,nd,additional_params)

if nargin<4
    
    additional_params=[];
    
end

% the parameters are ordered as
% [vec([C,B1,B2,...,Bp,SIG])',mc1_tp_i_j,...,mcn_tp_i_j]
% [vec([C,A0,A1,A2,...,Ap,S])',mc1_tp_i_j,...,mcn_tp_i_j]

% syntax for SVARs
% pnames=create_parameter_names(5,[0,3],1,{'chain',3},{'toto',5})

% syntax for VARs
% pnames=create_parameter_names(5,[1,3],1,{'chain',3},{'toto',5})

% improve the storing of variances and standard deviations

% assign each parameter to its governing chain and check that no parameter
% is assigned to multiple chains

% things are not ordered alphabetically as is done in RISE for all classes
% of models. If we decide to re-order things alphabetically...

is_svar=lags(1)==0;

if is_svar
    % svar
    prefix='a';
    
else
    % rfvar
    prefix='b';
    
end

p=vartools.param_creator(nvars,lags,constant,nd,additional_params,prefix,is_svar);

end