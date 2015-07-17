function [T1,W1]=newton_iteration_h(T0,dbf_plus,d0,dpb_minus,bf_cols,...
    pb_cols,kron_method,options)

% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

n=size(d0{1},1);
h=size(d0,2);
if nargin<7
    kron_method=[];
    if nargin<6
        pb_cols=[];
        if nargin<5
            bf_cols=[];
        end
    end
end

if isempty(kron_method)
    kron_method=true;
end
if isempty(bf_cols)
    bf_cols=1:n;
end
if isempty(pb_cols)
    pb_cols=1:n;
end
npb=numel(pb_cols);
if size(dbf_plus{1},2)~=numel(bf_cols)
    error('number of columns of dbf_plus inconsistent with the number of bf variables')
end
if size(dpb_minus{1},2)~=npb
    error('number of columns of dpb_minus inconsistent with the number of bp variables')
end

if isempty(T0)
    T0=zeros(n,npb,h);
end

n_npb=n*npb;
T0=reshape(T0,[n,npb,h]);
W=T0;
if kron_method
    G=zeros(n_npb*h);
else
    LMINUS=cell(1,h);
    LPLUS=cell(h);
end
% Lminus=zeros(npb);
Lplus01=zeros(n);
I_nx_nd=speye(n_npb);
for r0=1:h
    U=d0{r0};
    for r1=1:h
        U(:,pb_cols)=U(:,pb_cols)+dbf_plus{r0,r1}*T0(bf_cols,:,r1);
    end
    Ui=U\speye(n);
    T1_fi=-Ui*dpb_minus{r0};
    W(:,:,r0)=W(:,:,r0)-T1_fi;
    
    Lminus=-T1_fi(pb_cols,:);
    if kron_method
        Lminus_prime=Lminus.';
    else
        LMINUS{r0}=sparse(Lminus);
    end
    rows=(r0-1)*n_npb+1:r0*n_npb;
    for r1=1:h
        cols=(r1-1)*n_npb+1:r1*n_npb;
        Lplus01(:,bf_cols)=Ui*dbf_plus{r0,r1};
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
W=reshape(W,[n,npb*h]);
if kron_method
    % update T
    %---------
    delta=G\W(:);
else
%    delta=tfqmr(@(x)find_newton_step(x,LPLUS,LMINUS),-W(:),...
%        options.fix_point_TolFun);
	[delta,retcode]=transpose_free_quasi_minimum_residual(@(x)find_newton_step(x,LPLUS,LMINUS),-W(:),... 
			[],... %x0 initial guess
			options.fix_point_TolFun,... % tolerance level
			options.fix_point_maxiter,... % maximum number of iterations
			options.fix_point_verbose);
end

T1=T0+reshape(delta,[n,npb,h]);
if nargout>1
    W1=update_criterion();
end
T1=reshape(T1,[n,npb*h]);

    function W1=update_criterion()
        W1=T1;
        for r00=1:h
            U=d0{r00};
            for r11=1:h
                U(:,pb_cols)=U(:,pb_cols)+dbf_plus{r00,r11}*T1(bf_cols,:,r11);
            end
            Ui=U\eye(n);
            T1_fi=-Ui*dpb_minus{r00};
            W1(:,:,r00)=W1(:,:,r00)-T1_fi;
        end
    end

    function Gd=find_newton_step(delta,Lplus,Lminus)
        Gd=zeros(n*npb,h); % G*delta
        delta=reshape(delta,[n*npb,h]);
        for r00=1:h
            for r11=1:h
                Gd(:,r00)=Gd(:,r00)+vec(Lplus{r00,r11}*reshape(delta(:,r11),n,npb)*Lminus{r00});
            end
        end
        Gd=delta-Gd;
        Gd=Gd(:);
    end
end