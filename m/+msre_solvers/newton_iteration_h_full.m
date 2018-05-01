function [T1,W1]=newton_iteration_h_full(T0,Gplus01,A0,Aminus,kron_method,options)
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


[n,n,h]=size(A0);

if isempty(kron_method)
    kron_method=true;
end
if size(Gplus01,2)~=n
    error('number of columns of Gplus01 inconsistent with the number of variables')
end
if size(Aminus,2)~=n
    error('number of columns of Aminus inconsistent with the number of variables')
end

if isempty(T0)
    T0=zeros(n,n,h);
end

n_npb=n*n;
T0=reshape(T0,[n,n,h]);
W=T0;
if kron_method
    G=zeros(n_npb*h);
else
    LMINUS=cell(1,h);
    LPLUS=cell(h);
end
% Lminus=zeros(n);
% Lplus01=zeros(n);
I_nx_nd=speye(n_npb);
for r0=1:h
    U=A0(:,:,r0);
    for r1=1:h
        U=U+Gplus01(:,:,r0,r1)*T0(:,:,r1);
    end
    Ui=U\speye(n);
    T1_fi=-Ui*Aminus(:,:,r0);
    W(:,:,r0)=W(:,:,r0)-T1_fi;
    
    Lminus=-T1_fi;
    if kron_method
        Lminus_prime=Lminus.';
    else
        LMINUS{r0}=sparse(Lminus);
    end
    rows=(r0-1)*n_npb+1:r0*n_npb;
    for r1=1:h
        cols=(r1-1)*n_npb+1:r1*n_npb;
        Lplus01=Ui*Gplus01(:,:,r0,r1);
        if kron_method
            % build G
            %--------
            tmp=kron(Lminus_prime,Lplus01);
            if r0==r1
                tmp=tmp-I_nx_nd;
            end
            G(rows,cols)=tmp;
        else
            LPLUS{r0,r1}=sparse(Lplus01);
        end
    end
end
W=reshape(W,[n,n*h]);
if kron_method
    % update T
    %---------
    delta=G\W(:);
else
    delta0=[];
    [delta,retcode]=utils.optim.linear_systems_solver(...
        @(x)find_newton_step(x,LPLUS,LMINUS),-W(:),delta0,options);
end

T1=T0+reshape(delta,[n,n,h]);
if nargout>1
    W1=update_criterion();
end
T1=reshape(T1,[n,n*h]);

    function W1=update_criterion()
        W1=T1;
        for r00=1:h
            U=A0(:,:,r00);
            for r11=1:h
                U=U+Gplus01(:,:,r00,r11)*T1(:,:,r11);
            end
            Ui=U\eye(n);
            T1_fi=-Ui*Aminus(:,:,r00);
            W1(:,:,r00)=W1(:,:,r00)-T1_fi;
        end
    end

    function Gd=find_newton_step(delta,Lplus,Lminus)
        Gd=zeros(n*n,h); % G*delta
        delta=reshape(delta,[n*n,h]);
        for r00=1:h
            for r11=1:h
                Gd(:,r00)=Gd(:,r00)+vec(Lplus{r00,r11}*reshape(delta(:,r11),n,n)*Lminus{r00});
            end
        end
        Gd=delta-Gd;
        Gd=Gd(:);
    end
end