function kdata=embed(kdata,y,x)

nlags=kdata.nlags;

kdata.X=cell(nlags,1);

for ilag=1:nlags
    
    kdata.X{ilag}=y(:,(nlags+1:end)-ilag,:);
    
end

kdata.X=cell2mat(kdata.X);

if ~isempty(x)
    
    kdata.X=cat(1,x(:,nlags+1:end,:),kdata.X);
    
end

if kdata.nx*kdata.ng~=size(x,1)
    
    error('size problems in x probably due to initialization')
    
end

if kdata.nvars*kdata.ng~=size(y,1)
    
    error('mismatch between data and number of endogenous')
    
end

% kdata.nlags=nlags;

kdata.Y=y(:,nlags+1:end,:);

if any(isnan(kdata.X(:)))||any(isnan(kdata.Y(:)))
    
    error('nans in the data, you may want to change the estimation sample')
    
end

[kdata.K,kdata.T]=size(kdata.X);

end