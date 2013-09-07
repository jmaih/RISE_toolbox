function Aout=expand_array(Ain,flag)
h=numel(Ain);
[nrows,ncols]=size(Ain{1});
if flag==3
    h=ncols;
end
for st=1:h
    iter=(st-1)*nrows+1:st*nrows;
    switch flag
        case 1	% A0, Aminus or Aplus
            if st==1
                Aout=zeros(nrows*h);
            end
            Aout(iter,iter)=Ain{st};
        case 2 % B (shock impact)
            if st==1
                Aout=zeros(nrows*h,ncols);
            end
            Aout(iter,:)=Ain{st};
        case 3 % C (constant)
            if st==1
                Aout=zeros(nrows*ncols,1);
            end
            Aout(iter)=Ain{st};
    end
end
Aout={Aout};