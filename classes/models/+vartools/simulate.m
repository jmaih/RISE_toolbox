function [y,regimes]=simulate(y0,xdet,A,B,shocks,Qfunc,regimes)
% INTERNAL FUNCTION

if nargin<7

    regimes=[];

end

[nshocks,nperiods]=size(shocks);

if isempty(regimes)

    regimes=nan(1,nperiods);

end

if length(regimes)>nperiods

    regimes=regimes(1:nperiods);

elseif length(regimes)<nperiods

    miss=nperiods-length(regimes);

    regimes=[regimes(:).',nan(1,miss)];

end

nvars=nshocks;

y=zeros(nvars,nperiods);

y0t=y0;

ndet=size(xdet,2);

for t=1:nperiods

    [Q,retcode]=Qfunc(y0t(:,end));

    Q=full(Q);

    if t==1

        Qt=Q(:,:,ones(1,nperiods));

        h=size(A,3);

    else

        Qt(:,:,t)=Q;

    end

    rt=draw_regime();

    tdet=min(t,ndet);

    y(:,t)=A(:,:,rt)*flip_initial(y0t)+B(:,:,rt)*shocks(:,t);

    y0t=[y0t(:,2:end),y(:,t)];

end

    function rt=draw_regime()

        if isnan(regimes(t))

            if t==1
                % draw from initial distribution
                PAI=1/h*ones(1,h);

            else
                % draw conditional on yesterday's state
                PAI=Q(regimes(t-1),:);

            end

            cp=cumsum(PAI);

            cp=[0,cp(:).'];

            regimes(t)=find(cp>rand,1,'first')-1;

        end

        rt=regimes(t);

    end

    function xx=flip_initial(x)

        xx=fliplr(x); % xx=x(end:-1:1)

        xx=xx(:);

        if ~isempty(xdet)

            xx=[xdet(:,tdet);xx];

        end

    end

end