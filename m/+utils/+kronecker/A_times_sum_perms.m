function res=A_times_sum_perms(A,kron_matrices,matsizes,add_lead_term,varargin)
% A_times_sum_perms -- Computes the product between a matrix and a sum of
% tensor products
%
% ::
%
%
%   res=A_times_sum_perms(A,kron_matrices,matsizes,add_lead_term,P1,...,Pm)
%
% Args:
%
%    - **A** [matrix]:
%
%    - **kron_matrices** [cell array]: e.g. {M1,M2,...,Mn}, where the Mi's are
%    the matrices entering the various tensor products.
%
%    - **matsizes** [k x 2 vector]: sizes of the matrices. Note need to be
%    taken of the fact that these sizes are not necessarily the sizes of the
%    matrices entering **kron_matrices** . e.g. suppose we want to compute
%    A(UUB+UBU+BUU), where U is unknown but UU is known. Then we can have
%    kron_matrices={UU,B} but matsizes=[size(U);size(U);size(B)]
%
%    - **add_lead_term** [{true}|false]: if true, the term A(M1*M2*...*Mn) is
%    added to the sum
%
%    - **P1,...,Pm** [vectors]: permutations of the tensor products entering
%    the sum. The number of elements entering each Pi has to be the same as
%    the number of rows of **matsizes**, which represents the effective number
%    of matrices.
%
% Returns:
%    :
%
%    - **res** [matrix]: result of A(B*C*D+C*D*B+...)
%
% Note:
%
% Example:
%
%    See also:


nmat=size(matsizes,1); 

rc=prod(matsizes,1);

Ir=speye(rc(1));
Ic=speye(rc(2));

kron_order=fliplr(1:nmat);% reverse kronecker

rows_siz_flip=matsizes(kron_order,1).';
rows_span=reshape(1:rc(1),[1,rows_siz_flip]);

cols_siz_flip=matsizes(kron_order,2).';
cols_span=reshape(1:rc(2),[1,cols_siz_flip]);

A=sparse(A);
for imat=1:numel(kron_matrices)
        kron_matrices{imat}=sparse(kron_matrices{imat});
end
nterms=length(varargin);
algo=1;
switch algo 
    case 1
        % do the lead term
        %------------------
        if add_lead_term
            res=utils.kronecker.A_times_kron_Q1_Qk(A,kron_matrices{:});
        else
            res=sparse(size(A,1),rc(2));
        end
        % add the remaining ones
        %------------------------
        for ii=1:nterms
            res=res+do_one_permutation(varargin{ii});
        end 
    case 2
        ABC=utils.kronecker.kronall(kron_matrices{:});
        % do the lead term
        %------------------
        if add_lead_term
            res=A*ABC;
        else
            res=sparse(size(A,1),rc(2));
        end
        % add the remaining ones
        %------------------------
        for ii=1:nterms
            res=res+do_one_big_permutation(varargin{ii});
        end 
    case 3
        ABC=utils.kronecker.kronall(kron_matrices{:});
        
        % do the lead term
        %------------------
        if add_lead_term
            SumABC=ABC;
        else
            SumABC=sparse(size(ABC,1),size(ABC,2));
        end
        
        % add the remaining ones
        %------------------------
        for ii=1:nterms
            add_abc(varargin{ii});
        end 
        res=A*SumABC;
    case 4
        % vectorization strategy
        %------------------------
        ABC=utils.kronecker.kronall(kron_matrices{:});
        ncols=size(ABC,2);
        [nr,nc]=size(A);
        % do the lead term
        %------------------
        if add_lead_term
            SumR_ABC_L=kron(speye(ncols),A);
        else
            SumR_ABC_L=sparse(nr*ncols,nc*ncols);
        end
        
        % add the remaining ones
        %------------------------
        for ii=1:nterms
            [~,~,Li_,Ri_]=get_shufflers(varargin{ii});
            SumR_ABC_L=SumR_ABC_L+kron(Ri_.',A*Li_);
        end 
        res=reshape(SumR_ABC_L*ABC(:),nr,[]);
    otherwise
        error('unknown algorithm')
end

    function add_abc(in_order)
        newstyle=true;
        if newstyle
            [~,~,Li,Ri]=get_shufflers(in_order);
            SumABC=SumABC+Li*ABC*Ri;
        else
            [re_order_rows,re_order_cols]=get_shufflers(in_order);
            SumABC=SumABC+ABC(re_order_rows,re_order_cols);
        end
    end

    function [re_order_rows,re_order_cols,lfact,rfact]=get_shufflers(in_order)
        inflip=fliplr(in_order);
        new_order=kron_order(inflip);
        re_order_rows=reshape(permute(rows_span,[1,new_order+1]),1,[]); 
        re_order_cols=reshape(permute(cols_span,[1,new_order+1]),1,[]); 
        % first step: permute rows
        %-------------------------
        % Note that A*factor is different from A(:,re_order_rows) !!!
        if nargout>2
            lfact=Ir(re_order_rows,:);
            if nargout>3
                rfact=Ic(:,re_order_cols);
            end
        end
    end

    function res=do_one_big_permutation(in_order)
        [re_order_rows,re_order_cols]=get_shufflers(in_order);  
        res=A*ABC(re_order_rows,re_order_cols);
    end

    function res=do_one_permutation(in_order)
        [~,re_order_cols,lfact]=get_shufflers(in_order);  
        % first step: permute rows
        %-------------------------
        % Note that A*lfact is different from A(:,re_order_rows) !!!
        res=utils.kronecker.A_times_kron_Q1_Qk(A*lfact,kron_matrices{:});

        % second step: permute columns
        %-----------------------------
        % the two options below seem to take about the same time to
        % compute...
        res=res(:,re_order_cols);% res=res*rfact;
        
%         ABC=utils.kronecker.kronall(kron_matrices{:});
%         Ic=speye(rc(2));
%         reo=kron_matrices(in_order);
%         ACB=utils.kronecker.kronall(reo{:});
%         max(max(abs(ACB-ABC(re_order_rows,re_order_cols))))
%         max(max(abs(ACB-Ir(re_order_rows,:)*ABC*Ic(:,re_order_cols))))
%         max(max(abs(A*ACB-res)))
    end
end

%{
clc
kronall=@utils.kronecker.kronall;
imax=10;
imin=3;
matrices={
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    };
nmat=numel(matrices);

matsizes=zeros(nmat,2);

for ii=1:nmat
    matsizes(ii,:)=size(matrices{ii});
end

A=rand(20,prod(matsizes(:,1)));

orders={randperm(nmat),randperm(nmat),randperm(nmat),randperm(nmat),randperm(nmat)};

no=numel(orders);

tic
Fasit=A*kronall(matrices{:});
for ialt=1:no
    tmp=matrices(orders{ialt});
    Fasit=Fasit+A*kronall(tmp{:});
end
toc

tic
res=utils.kronecker.A_times_sum_perms(A,matrices,matsizes,true,...
    orders{:});
toc

max(abs(Fasit(:)-res(:)))

%}