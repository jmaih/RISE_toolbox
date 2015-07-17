function oo = dv_vz_omega(dv_vz,nz,index)

nd=size(dv_vz,1);

j =  2;
k =  3;
l =  4;
m =  5;
n =  6;
switch index
    case 1
        layers=3;
        oo=isum([k,l,j],[j,l,k],[j,k,l]);
    case 2
        layers=4;
        oo=isum([l,m,k,j],[k,m,l,j],[j,m,l,k],[k,l,m,j],...
            [j,l,m,k],[j,k,m,l]);
    case 3
        layers=4;
        oo=isum([k,l,m,j],[j,l,m,k],[j,k,m,l],[j,k,l,m]);
    case 4
        layers=4;
        oo=isum([k,l,j,m],[j,l,k,m],[j,k,l,m]);
    case 5
        layers=5;
        oo=isum([m,n,l,k,j],[l,n,m,k,j],[k,n,m,l,j],...
            [j,n,m,l,k],[l,m,n,k,j],[k,m,n,l,j],...
            [j,m,n,l,k],[k,l,n,m,j],[j,l,n,m,k],...
            [j,k,n,m,l]);
    case 6
        layers=5;
        oo=isum([l,m,n,k,j],[k,m,n,l,j],[j,m,n,l,k],...
            [k,l,n,m,j],[j,l,n,m,k],[j,k,n,m,l],...
            [k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [j,k,l,n,m]);
    case 7
        layers=5;
        oo=isum([l,m,k,n,j],[k,m,l,n,j],[j,m,l,n,k],...
            [k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [l,m,j,n,k],[k,m,j,n,l],[j,m,k,n,l],...
            [k,l,j,n,m],[j,l,k,n,m],[j,k,l,n,m],...
            [k,l,j,m,n],[j,l,k,m,n],[j,k,l,m,n]);
    case 8
        layers=5;
        oo=isum([k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [j,k,l,n,m],[j,k,l,m,n]);
    case 9
        layers=5;
        oo=isum([k,l,m,j,n],[j,l,m,k,n],[j,k,m,l,n],...
            [j,k,l,m,n],[k,l,n,j,m],[j,l,n,k,m],...
            [j,k,n,l,m],[j,m,n,k,l],[k,m,n,j,l],[l,m,n,j,k]);
    otherwise
        error('index must be in [1,9]')
end

    function ss=isum(varargin)
        if issparse(dv_vz)
            dv_vz=full(dv_vz);
        end
        dv_vz=reshape(dv_vz,[nd,nz*ones(1,layers)]);
        ss=0;
        for ii=1:length(varargin)
            ss=ss+ipermute(dv_vz,[1,varargin{ii}]);
        end
        ss=reshape(ss,[nd,nz^layers]);
        ss=sparse(ss);
    end
end

