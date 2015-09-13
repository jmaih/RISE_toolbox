function [dd,di]=tnrnd(a,b,m,s,options)
if nargin<5
    options=[];
    if nargin<4
        s=[];
        if nargin<3
            m=[];
        end
    end
end

if isempty(options)
    options=utils.forecast.rscond.tmvnrnd();
end

if isempty(s)
    s=1;
end
if ~isfinite(s)||~isreal(s)||s<0
    error('standard deviation must be real and >=0')
end
if isempty(m)
    m=0;
end

r=rand;

dd=direct_draw();
if nargout>1
    di=indirect_draw();
end
    function d=direct_draw()
        d=m;
        if s>=0
            Fa=normcdf(a,m,s);
            Fb=normcdf(b,m,s);
            u=Fa+r*(Fb-Fa);
            d=apply_correction(norminv(u,m,s));
        end
    end

    function d=indirect_draw()
        d=m;
        if s>=0
            ai=(a-m)/s;
            bi=(b-m)/s;
            Fa=normcdf(ai);
            Fb=normcdf(bi);
            u=Fa+r*(Fb-Fa);
            d=norminv(u);
            d=apply_correction(m+s*d);
        end
    end
    function d=apply_correction(d0)
        d=max(d0,a);
        d=min(d,b);
        if options.debug && d~=d0
            disp([d,d0])
        end
    end
end