function [yk,betta,pcfit]=three_pass_regression_filter(y,X,Z)
% inputs: 
%       y is a T x 1 vector to be predicted
%       X is a T x N matrix of predictors
%       Z is a T x M matrix of proxies. If Z is a scalar, then M is set to Z and
% the proxies are constructed automagically.
% y, X, and Z could be rise_time_series objects
% outputs:
%         yk: predicted values
%         betta: I forgot, see paper
%         pcfit: principal component regression fit

% k=4;
% figure;
% for ii=1:k
%     [yk,betta,pcr_fit]=three_pass_regression_filter(yy,xx,2+ii-1);
%     subplot(2,2,ii)
%     plot(m_gdp+[yy(300:end),yk(300:end),pcr_fit(300:end)]),legend({'actual',['3prf(',int2str(2+ii-1),')'],'pcr'})
% end
if nargin<3
    Z=[];
end
is_time_series=strcmp(class(y),'rise_time_series');
if is_time_series && ~strcmp(class(X),'rise_time_series')
    error([mfilename,':: when one input is a time series, so must all inputs'])
end
automatic=isempty(Z)||isscalar(Z);
if is_time_series && ~isempty(Z) && ~isscalar(Z)
    if ~strcmp(class(Z),'rise_time_series')
        error([mfilename,':: when one input is a time series, so must all inputs'])
    end
end    
if is_time_series 
    last_date=y.finish;
    [y,X]=intersect(y,X);
    if ~automatic
        [y,Z]=intersect(y,Z);
        [X,Z]=intersect(Z,Z);
        Z=double(Z);
    end
    TimeInfo=y.TimeInfo;
    if ~strcmp(last_date,y.finish)
        error([mfilename,':: last date is lost during the intersection procedure'])
    end
    y=double(y);
    X=double(X);
end
if isnan(y(end))||any(isnan(X(end,:)))||(~automatic && any(isnan(Z(end,:))))
    error([mfilename,':: for forecasting purposes, last observation cannot contain nans'])
end

[T,N]=size(X);
if size(y,1)~=T
    error([mfilename,':: matrix dimensions are inconsistent'])
end
%% Finding proxies: Statistical proxy-Selection algorithm
if automatic
    algo=any(isnan(y))||any(any(isnan(X)));
    % set M
    if isscalar(Z)
        M=Z;
    else
        M=10;
    end
    
    r0=y;
    Z=nan(T,M);
    for k=1:M
        Z(:,k)=r0;
        if algo==0
            if k==1
                JT=eye(T)-1/T*ones(T);
                JN=eye(N)-1/N*ones(N);
                ybar=1/T*ones(1,T)*y;
            end
            [yk,betta]=three_pass_no_nans(Z(:,1:k));
        else
            [yk,betta]=three_pass_with_nans(Z(:,1:k));
        end
        r0=y-yk;
    end
else
    algo=any(isnan(y))||any(any(isnan(X)))||any(any(isnan(Z)));
    if algo==0
        JT=eye(T)-1/T*ones(T);
        JN=eye(N)-1/N*ones(N);
        ybar=1/T*ones(1,T)*y;
        [yk,betta]=three_pass_no_nans(Z);
    else
        [yk,betta]=three_pass_with_nans(Z);
    end
end
if is_time_series 
    yk=rise_time_series(TimeInfo,yk);
end
if nargout>2
    pcfit=pcr();
    if is_time_series
        pcfit=rise_time_series(TimeInfo,pcfit);
    end
end

    function [yhat,betta]=three_pass_no_nans(Z)    
        alpha=JN*X'*JT*Z/...
            (Z'*JT*X*JN*X'*JT*X*JN*X'*JT*Z)...
            *Z'*JT*X*JN*X'*JT*y;
        yhat=ybar+JT*X*alpha;
        if nargout>1
            betta=(Z'*JT*Z)\Z'*JT*X*alpha;
        end
    end

    function [yhat,betta]=three_pass_with_nans(Z)
        % first pass
        PHI=nan(M+1,N);
        rhs=[ones(T,1),Z];
        for ii=1:N
            lhs=X(:,ii);
            PHI(:,ii)=ols(lhs,rhs);
        end
        % second pass
        rhs=[ones(N,1),transpose(PHI(2:end,:))];
        F=nan(M+1,T);
        for t=1:T
            lhs=transpose(X(t,:));
            F(:,t)=ols(lhs,rhs);
        end
        % third pass
        rhs=[ones(T,1),transpose(F(2:end,:))];
        lhs=y;
        betta=ols(lhs,rhs);
        yhat=rhs*betta;
    end

    function pcfit=pcr()
        % Principal component Analysis
        good=~isnan(y) & ~any(isnan(X),2);
        
        normalize=inline('bsxfun(@rdivide,bsxfun(@minus,x,mean(x)),std(x))','x');
        Xbar=normalize(X(good,:));
        [V,D] = eig(Xbar'*Xbar);
        D=diag(D);
        [~,tag]=sort(abs(D),1,'descend');
        D=D(tag);
        V=V(:,tag);
        XX0=0;
        for ii=1:numel(D)
            XX0=XX0+D(ii)*(V(:,ii)*V(:,ii)');
        end
%         % test that we get back what we should...
%         max(max(abs(XX0-Xbar'*Xbar)))
        % Now we construct the factors
        F=nan(T,M);
        for m=1:M
            F(:,m)=X*V(:,m);
        end
        % now we can regress
        pcfit=F*(F\y);
    end
end

function b=ols(y,X)
good=~isnan(y);
good=good & ~any(isnan(X),2);
b=X(good,:)\y(good);
end

%{


plot((1:T)',[y,yhat])
%% The two procedures produce the same answer for a given Z
max(abs(yhat-yk))

%% Principal component Analysis
normalize=inline('bsxfun(@rdivide,bsxfun(@minus,x,mean(x)),std(x))','x');
Xbar=normalize(X);
[V,D] = eig(Xbar'*Xbar);
D=diag(D);
[junk,tag]=sort(abs(D),1,'descend');
D=D(tag);
V=V(:,tag);
XX0=0;
for ii=1:numel(D)
    XX0=XX0+D(ii)*(V(:,ii)*V(:,ii)');
end
max(max(abs(XX0-Xbar'*Xbar)))
% Now we construct the factors
F=nan(T,M);
for m=1:M
    F(:,m)=X*V(:,m);
end
% now we can regress
pcfit=F*(F\y);
%% plot all
plot((1:T)',[y,yk,yhat,pcfit])
legend('y','3prf1','3prf2','pcfit')
%}
