function [p,priors]=create_parameters_switch(start_from_mode)

[p,priors]=create_parameters(start_from_mode);

p=rmfield(p,'sig_Eps_phi');

priors=rmfield(priors,'sig_Eps_phi');

priors.sig_Eps_phi_hetero_1 ={ 0.01, 0.0001, 2, 'inv_gamma_pdf(0.9)', 0.000000000001, 100};

priors.sig_Eps_phi_hetero_2 ={ 0.01, 0.0001, 2, 'inv_gamma_pdf(0.9)', 0.000000000001, 100};

priors.hetero_tp_1_2    ={ 0.7  , 0.0256, 0.7761 , 'beta_pdf(0.9)'  , .00000001, .99999999999};

priors.hetero_tp_2_1    ={ 0.7  , 0.0256, 0.7761 , 'beta_pdf(0.9)'  , .00000001, .99999999999};

if start_from_mode
    
    priors.sig_Eps_phi_hetero_1{1}=   0.03;
    
    priors.sig_Eps_phi_hetero_2{1}=   0.08;
    
    priors.hetero_tp_1_2{1} =1-0.9794;
    
    priors.hetero_tp_2_1{1} =1-0.9662;

end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);

for ip=1:numel(fields)
    
    name=fields{ip};
    
    p.(name)=priors.(name){1};
    
end

end