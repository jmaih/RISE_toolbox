function oo_=moments(db,drange,ar,lambda)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<4
    
    lambda=[];
    
    if nargin<3
        
        ar=[];
        
        if nargin<2
            
            drange=[];
            
        end
        
    end
    
end

if isempty(ar)
    
    ar=1;
    
end

if ~isempty(drange)
    
    db=db(drange);
    
end

oo_.mean=mean(db);

if ~isempty(lambda)
    
    [oo_.hp_trend,oo_.hp_cycle] = hpfilter(db,lambda);
    
    db=oo_.hp_cycle;
    
else
    
    db = bsxfun(db,@minus,oo_.mean);
    
end

oo_.vcov=cov(db);

oo_.skewness=skewness(db);

oo_.kurtosis=kurtosis(db);

oo_.variance = diag(oo_.vcov);oo_.variance=oo_.variance(:).';

oo_.stdev = sqrt(oo_.variance);

oo_.corrcoef = corrcoef(db);

if ar > 0
    
    y=double(db);
    
    for i=1:ar
        
        oo_.autocorr{i} = y(ar+1:end,:)'*y(ar+1-i:end-i,:)./((size(y,1)-ar)*std(y(ar+1:end,:))'*std(y(ar+1-i:end-i,:)));
        
    end
    
end

end