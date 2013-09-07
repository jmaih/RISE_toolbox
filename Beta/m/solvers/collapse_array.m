function Aout=collapse_array(Ain,flag,hh)
[nrows,ncols,order]=size(Ain);
nvar=nrows/hh;

for st=1:hh
    iter=(st-1)*nvar+1:st*nvar;
    switch flag
        case 1
            if st==1
                Aout=zeros(nvar,nvar,hh);
            end
            Aout(:,:,st)=Ain(iter,iter);
        case 2
            if st==1
                Aout=zeros(nvar,ncols,order,hh);
            end
			for oo=1:order
	            Aout(:,:,oo,st)=Ain(iter,:,oo);
			end
        case 3
            if st==1
                Aout=zeros(nvar,hh);
            end
            Aout(:,st)=Ain(iter);
    end
end
