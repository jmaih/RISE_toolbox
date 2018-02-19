function kdata=embed(kdata,y,x)

nlags=kdata.nlags;

kdata.estim_.X=cell(nlags,1);

for ilag=1:nlags
    
    kdata.estim_.X{ilag}=y(:,(nlags+1:end)-ilag,:);
    
end

kdata.estim_.X=cell2mat(kdata.estim_.X);

if ~isempty(x)
    
    kdata.estim_.X=cat(1,x(:,nlags+1:end,:),kdata.estim_.X);
    
end

if kdata.nx*kdata.ng~=size(x,1)
    
    error('size problems in x probably due to initialization')
    
end

if kdata.nvars*kdata.ng~=size(y,1)
    
    error('mismatch between data and number of endogenous')
    
end

% kdata.nlags=nlags;

kdata.estim_.Y=y(:,nlags+1:end,:);

if any(isnan(kdata.estim_.X(:)))||any(isnan(kdata.estim_.Y(:)))
    
    error('nans in the data, you may want to change the estimation sample')
    
end

[kdata.estim_.K,kdata.estim_.T]=size(kdata.estim_.X);

end