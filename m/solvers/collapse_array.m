function Aout=collapse_array(Ain,flag,hh)
% takes as input a multidimensional array and outputs a cell array
[nrows,ncols,order]=size(Ain);
nvar=nrows/hh;

Aout=cell(1,hh);
for st=1:hh
    iter=(st-1)*nvar+1:st*nvar;
    switch flag
        case 1
            Aout{st}=Ain(iter,iter);
        case 2
            Aout{st}=zeros(nvar,ncols,order);
            for oo=1:order
                Aout{st}(:,:,oo)=Ain(iter,:,oo);
            end
        case 3
            Aout{st}=Ain(iter);
    end
end
