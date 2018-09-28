function names=perturbation_control_param_names(parameterize)
% name for the perturbation parameter entering the structural model
% explicitly
% if parameterize is true, return also the default (MAIH's perturbation)
% values

if nargin==0
    
    parameterize=false;
    
end

tail='___';

names=strcat({'sig','iota_1','iota_2'},tail);

if ~parameterize
    
    return
    
end

sig_val=0;

names=[names
    {sig_val,1,0}];

names=names.';   

end