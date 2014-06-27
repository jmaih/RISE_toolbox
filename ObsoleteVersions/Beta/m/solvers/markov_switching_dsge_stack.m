function [A0,C,Aminus,Aplus,B]=markov_switching_dsge_stack(Aminus_st,A0_st,Aplus_st,B_st,...
    C_st,... % endo_nbr x h matrix of constant
    P,... % h x h transition matrix (Pij=prob(i-->j)) such that sum_j(Pij)=1
    T,... % endo_nbr x endo_nbr x h array of solution..
    delta_st... % endo_nbr x h matrix of solution constants
    )
if nargin<4||nargin>8
    error([mfilename,':: wrong number of arguments'])
end

[endo_nbr,exo_nbr,h]=size(B_st);

if nargin<8
    delta_st=[];
    if nargin<7
        T=[];
        if nargin<6
            P=1;
            if nargin<5
                C_st=[];
            end
        end
    end
end
nstar=endo_nbr*h;

Args={'Aminus_st','A0_st','Aplus_st','B_st','C_st','P','T','delta_st'
    [endo_nbr,endo_nbr,h],[endo_nbr,endo_nbr,h],[endo_nbr,endo_nbr,h],...
    [endo_nbr,exo_nbr,h],[endo_nbr,h],[h,h],[nstar,nstar],[nstar,1]};
for ii=1:size(Args,2)
    eval([Args{1,ii},'=CheckArgument(',Args{1,ii},',Args{2,ii},ii);'])
end

if any(sum(P,2)~=1)
    Ptmp=transpose(P);
    if all(sum(Ptmp,2)==1)
        P=Ptmp;
    else
        error([mfilename,':: columns of transition matrix should sum to 1'])
    end
end

if nargout>0
    Aminus=zeros(nstar);
    C=zeros(nstar,1);
    if nargout>2
        A0=zeros(nstar);
        Aplus=zeros(nstar);
        B=zeros(nstar,exo_nbr);
    end
end

% this seems faster than calling setdfiff
regimes=1:h;
others=cell(1,h);
for st=1:h
    others{st}=[regimes(1:st-1),regimes(st+1:end)];
end
for st=1:h
    iter=(st-1)*endo_nbr+1:st*endo_nbr;
    if nargout>0
        A0(iter,iter)=A0_st(:,:,st);
        C(iter,1)=C_st(:,st);
        if h>1
            tmp=0;
            tmp_c=0;
            for j=others{st} % setdiff(1:h,st)
                it_=(j-1)*endo_nbr+1:j*endo_nbr;
                tmp=tmp+P(st,j)*T(it_,it_);
                tmp_c=tmp_c+P(st,j)*delta_st(it_,1);
            end
            A0(iter,iter)=A0(iter,iter)+Aplus_st(:,:,st)*tmp;
            C(iter,1)=C(iter,1)+Aplus_st(:,:,st)*tmp_c;
        end
        if nargout>2
            B(iter,:)=B_st(:,:,st);
            Aplus(iter,iter)=P(st,st)*Aplus_st(:,:,st);
            Aminus(iter,iter)=Aminus_st(:,:,st);
        end
    end
end

