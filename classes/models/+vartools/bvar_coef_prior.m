function B=bvar_coef_prior(coefprior,nvars,nlags,nx)

if isempty(coefprior)
    
    coefprior=eye(nvars);
    
elseif isscalar(coefprior)
    
    coefprior=diag(coefprior*ones(1,nvars));
    
elseif numel(coefprior)==nvars
    
    coefprior=diag(coefprior);
    
end

if ~isequal(size(coefprior),[nvars,nvars])
    
    error('wrong specification of coefprior')
    
end


B=zeros(nvars,nx);
    

B=[B,coefprior,zeros(nvars,(nlags-1)*nvars)];

end