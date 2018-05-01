function [iterating_function,...
solution_function,...
inverse_solution_function,...
sampling_function] = fwz_newton_system(Gplus01,A0,Aminus,Q)
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

% We need to rewrite the system in the form
% A(st)*x=B(st)*x_{t-1}+PSI*e_t+PI(st), with the ll forward-looking
% variables at the bottom
% just conforming to the FWZ notation. In this case, the columns will sum
% to 1.
% N.B: The case where A, as constructed below, is singular is not
% implemented.
hh=size(Q,1);
Aplus=unearth_frwrd_matrix();
 
fwl=any(Aplus(:,:,1));
for istate=2:hh
    fwl=fwl|any(Aplus(:,:,istate));
end

ll=sum(fwl);
nn=size(Aplus,1)+ll;
Inl=eye(nn-ll);
Il=eye(ll);
Inl0=[Inl,zeros(nn-ll,ll)];
ZnlIl=[zeros(nn-ll,ll);-Il];

A=zeros(nn,nn,hh);
B=zeros(nn,nn,hh);
PIE=[zeros(nn-ll,ll);Il]; % SAME common across states
for istate=1:hh
    A(1:nn-ll,:,istate)=[A0(:,:,istate),Aplus(:,fwl,istate)];
    A(nn-ll+(1:ll),fwl,istate)=eye(ll);
    B(1:nn-ll,1:nn-ll,istate)=-Aminus(:,:,istate);
    B(nn-ll+(1:ll),nn-ll+(1:ll),istate)=eye(ll);
end

solution_function=@final_solution;
iterating_function=@X_update;
inverse_solution_function=@fwz_T_into_X;
sampling_function=@sample_guess;

    function [X,T]=sample_guess()
        X=zeros(ll,nn-ll,hh);
        for ii=1:hh
            Xbar=randn(nn);
            [QQ,RR]=qr(Xbar); %#ok<NASGU>
            V=QQ(:,1:nn-ll);
            AV=A(:,:,ii)*V;
            IX=AV;
            X(:,:,hh)=-IX(nn-ll+1:end,:);
%             I do not understand the formula below
%             X(:,:,hh)= V*inv(AV(1:nn-ll,1:nn-ll));
        end
        T=final_solution(X);
    end

    function X=fwz_T_into_X(T)
        X=zeros(ll,nn-ll,hh);
        for s_now=1:hh
            T2=0;
            for s_lead=1:hh
                T2=Q(s_now,s_lead)*T(:,:,s_lead);
            end
            T2=T2*T(:,:,s_now);
            Gi=[T(:,:,s_now);T2(fwl,:)];
            [Qi,Ri]=qr(Gi); %#ok<NASGU>
            Vi=Qi(:,1:nn-ll);
            IX=A(:,:,s_now)*Vi;
            X(:,:,s_now)=-IX(nn-ll+1:end,:);
        end
    end

    function [X,F]=X_update(X)
        if isempty(X)
            X=zeros(ll,nn-ll,hh);
        elseif isequal(size(X),[nn-ll,nn-ll,hh])
            X=fwz_T_into_X(X);
        end
        FPRIME=zeros(hh*ll*(nn-ll));
        F=zeros(ll,hh*(nn-ll));
        for s_now=1:hh
            XBA=0;
            now_rows=(s_now-1)*ll*(nn-ll)+1:s_now*ll*(nn-ll);
            f_rows=(s_now-1)*(nn-ll)+1:s_now*(nn-ll);
            for s_plus=1:hh
                plus_cols=(s_plus-1)*ll*(nn-ll)+1:s_plus*ll*(nn-ll);
                
                BoA=B(:,:,s_plus)/A(:,:,s_now);
                
                FPRIME(now_rows,plus_cols)=Q(s_now,s_plus)*kron(transpose(Inl0*BoA*[Inl;-X(:,:,s_now)]),Il);
                
                XI=[X(:,:,s_plus),Il];
                
                XBA=XBA+Q(s_now,s_plus)*XI*BoA*ZnlIl;
                
                F(:,f_rows)=F(:,f_rows)+Q(s_now,s_plus)*XI*BoA*[Inl;-X(:,:,s_now)];
            end
            FPRIME(now_rows,now_rows)=FPRIME(now_rows,now_rows)+kron(Inl,XBA);
        end
        X=reshape(X,ll,(nn-ll)*hh);
        X=X(:)-FPRIME\F(:);
        X=reshape(X,[ll,nn-ll,hh]);
    end

    function G=X2Gamma(X)
        G=zeros(nn,nn,hh);
        for kk=1:hh
            F12=[[Inl;-X(:,:,kk)],PIE]\B(:,:,kk); % from equation (7)
            V=A(:,:,kk)\[Inl;-X(:,:,kk)]; % from equation (11)
            G(:,:,kk)=V*F12(1:nn-ll,:); % from equation (5)
        end
    end

    function T=final_solution(X)
        GAM=X2Gamma(X);
        T=GAM(1:nn-ll,1:nn-ll,:);
    end

    function Aplus=unearth_frwrd_matrix()
        % this function extracts Aplus from Gplus and is used for solving models
        % using the Waggoner-Zha approach, which states the problem to solve as
        % Aplus(st)*X_{t+1}+A0(st)*X_{t}+Aminus(st)*X_{t-1}+B(st)*Et=0 whereas in
        % my approach, the problem solved would be
        % Aplus(st+1)*X_{t+1}+A0(st)*X_{t}+Aminus(st)*X_{t-1}+B(st)*Et=0
        endo_nbr=size(Gplus01,1);
        Aplus=zeros(endo_nbr,endo_nbr,hh);
        for reg=1:hh
            Aplus(:,:,reg)=Gplus01(:,:,reg,reg)/Q(reg,reg);
        end
    end
end

