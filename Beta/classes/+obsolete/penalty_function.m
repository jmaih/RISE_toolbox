function p=penalty_function(gx,k,type)
alpha_=100;
d=100;
if nargin<3
    type='';
    if nargin<2
        k=[];
    end
end

if isempty(type)
    if isempty(k)
        type='static';
    else
        type='dynamic';
    end
end

gx=gx(:);
xi=max(0,gx);
switch type
    case 'death'
        if any(xi)
            p=inf;
        end
    case 'static'
        p=d*sum(xi);
    case 'dynamic'
        dt=k^.1;
        p=dt*sum(theta_xi().*xi.^gamma_xi());
end

    function g=gamma_xi
        g=ones(size(xi));
        g(xi>1)=2;
    end

    function thet=theta_xi
        thet=zeros(size(xi));
        thet(xi<1e-5)=alpha_;
        thet(xi>=1e-5 & xi<1e-3)=10*alpha_;
        thet(xi>=1e-3 & xi<1)=100*alpha_;
        thet(xi>1)=1000*alpha_;
    end
end