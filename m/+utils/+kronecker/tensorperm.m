function [inew,jnew]=tensorperm(matsizes,varargin)
% tensorperm -- finds indices (rows and columns) for the permutations of a
% tensor
%
% ::
%
%
%   [irows,jcols]=tensorperm(matsizes,order_1,order_2,...,order_n)
%
%   [irows,jcols]=tensorperm(matsizes,order_1,order_2,...,order_n,'grid')
%
% Args:
%
%    - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
%    entering the tensor. Each row represents the size of a matrix and it is
%    assumed that the main(or first) tensor product is ordered [1,2,...,k]
%
%    - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
%    1,2,...,k should be part of the vector and the vector should have exactly
%    k elements
%
%    - **grid** [string]: if present, a grid is used to compute the indexes of
%    the main kronecker product
%
% Returns:
%    :
%
%    - **irows** [vector|matrix]: each column corresponds to the rows for a
%    particular ordering as listed in varargin
%
%    - **jcols** [vector|matrix]: each column corresponds to the columns for a
%    particular ordering as listed in varargin
%
% Note:
%
% Example:
%
%    See also: sum_permutations

nargs=length(varargin);

grid_it=false;

if ischar(varargin{end})
    
    grid_it=true;
    
    varargin=varargin(1:end-1);

    nargs=nargs-1;

end

nmat=size(matsizes,1);

inew=do_cols(1);

jnew=do_cols(2);

    
    function index=do_cols(ref)
        
        siz=matsizes(:,ref).';
        
        nrows=prod(siz);
        
        if grid_it
        
            old_=utils.gridfuncs.mygrid(siz);
        
        else
            
            new_order=1:nmat;
        
            old_=find_indexes();
        
        end
        
        index=zeros(prod(siz),nargs);
        
        for iarg=1:nargs
            
            new_order=varargin{iarg};
            
            coefs=compute_factors();
            
            tmp=old_(:,new_order(end));% new_old_=old_(:,new_order); tmp=new_old_(:,end);
            
            for icol=1:nmat-1
            
                tmp=tmp+(old_(:,new_order(icol))-1)*coefs(icol); %tmp=tmp+(new_old_(:,icol)-1)*coefs(icol);
            
            end
            
            index(tmp,iarg)=1:numel(tmp);
        
        end
        
        function fact=compute_factors()
            
            fact=siz(new_order);
            
            fact=fliplr([fact(2:end),1]);
        
            fact=fliplr(cumprod(fact));
        
        end
        
        function ind=find_indexes()
            
            rowcol=compute_factors();
            
            ind=nan(nrows,nmat);
            
            x=(1:nrows).';
            
            for icol_=1:nmat-1
            
                ind(:,icol_)=iterate(rowcol(icol_));
            
            end
            
            ind(:,end)=x;
            
            function ind=iterate(y)
                
                [rest,n]=mymod(x-1,y); % <--rest=mod(x-1,y);  n=(x-1-rest)/y;
                
                ind=n+1;
            
                x=rest+1;
        
            end
            
        end
        
    end

end

function [rest,n]=mymod(x,y)

n=floor(x./y);

rest=x-y.*n;

end