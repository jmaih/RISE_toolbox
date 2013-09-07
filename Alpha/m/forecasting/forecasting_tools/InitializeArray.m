function Out=InitializeArray(fcn,nrows,ncols,npages,x0)
if nargin<5
    x0=[];
    if nargin<4
        npages=1;
        if nargin<3
            error([mfilename,':: Too few arguments'])
        end
    end
end
nlags=size(x0,2);
if isempty(x0)
    x0=zeros(nrows,nlags);
end
Out=fcn(nrows,ncols+nlags,npages);
Out(:,1:nlags,:)=x0(:,:,ones(npages,1));
end
