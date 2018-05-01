function db=prctile(db,p)
% prctile Percentiles of a time series (ts)
%
% ::
%
%
%   db=prctile(db,p)
%
% Args:
%
%    - **db** [ts] : time series with many pages (third dimension). The time
%      series may have one or several variables.
%
%    - **p** [scalar|vector] : scalar or a vector of percent values
%
% Returns:
%    :
%
%    - **db** [ts] : time series with as many pages as the length of **p**.
%
% Note:
%
% Example:
%
%    test=ts(1990,rand(100,3,200),{'a','b','c'});
%    tmp=prctile(test,[10,50,90])
%    plot(tmp('a'))
%
%    See also:

if ~isvector(p) || numel(p) == 0 || any(p < 0 | p > 100) || ~isreal(p)
    
    error('percentiles must be between 0 and 100');
    
end

x=double(db);

if size(x,3)==1
    
    x=permute(x,[1,3,2]);
    
    if ~isempty(db.varnames{1})
        
        error('names found in the database but only one page. Suppress names or add pages')
        
    end
    
end

[nobs,nvars,npages]=size(x);

np=numel(p);

y=nan(nobs,nvars,np);

for ivar=1:nvars
    
    y(:,ivar,:) = prctile_engine(permute(x(:,ivar,:),[1,3,2]));
    
end

if isempty(db.varnames)
    
    db=ts(db.start,y);
    
else
    
    db=ts(db.start,y,db.varnames);
    
end

    function y = prctile_engine(x)
        
        perm = [2 1];
        x = permute(x,perm);
        
        x = sort(x,1);
        nonnans = ~isnan(x);
        
        % If there are no NaNs, do all cols at once.
        if all(nonnans(:))
            n = npages;
            if isequal(p,50) % make the median fast
                if rem(n,2) % n is odd
                    y = x((n+1)/2,:);
                else        % n is even
                    y = (x(n/2,:) + x(n/2+1,:))/2;
                end
            else
                q = [0 100*(0.5:(n-0.5))./n 100]';
                xx = [x(1,:); x(1:n,:); x(n,:)];
                y = zeros(numel(p), nobs);
                y(:,:) = interp1q(q,xx,p(:));
            end
            
            % If there are NaNs, work on each column separately.
        else
            % Get percentiles of the non-NaN values in each column.
            y = nan(numel(p), nobs);
            for j = 1:nobs
                nj = find(nonnans(:,j),1,'last');
                if nj > 0
                    if isequal(p,50) % make the median fast
                        if rem(nj,2) % nj is odd
                            y(:,j) = x((nj+1)/2,j);
                        else         % nj is even
                            y(:,j) = (x(nj/2,j) + x(nj/2+1,j))/2;
                        end
                    else
                        q = [0 100*(0.5:(nj-0.5))./nj 100]';
                        xx = [x(1,j); x(1:nj,j); x(nj,j)];
                        y(:,j) = interp1q(q,xx,p(:));
                    end
                end
            end
        end
        
        % undo the DIM permutation
        y = ipermute(y,perm);
        
    end
end