function [tr,utr]=bounds_transform(lb,ub)
%
% transforms constrained parameters into unconstrained ones for use in root
% finders or optimization routines that do not handle constraints on
% parameters.
%
% returns function handles
% - tr: transform parameters into unconstrained units
% - utr: unstransform parameters into original units
%
% ------------------------------- tests ----------------------------------
% [tr,utr]=utils.estim.bounds_transform(0,inf),abs(utr(tr(11))-11)<1e-13
% [tr,utr]=utils.estim.bounds_transform(-inf,0),abs(utr(tr(-11))+11)<1e-13
% [tr,utr]=utils.estim.bounds_transform(10,15),abs(utr(tr(11))-11)<1e-13
% [tr,utr]=utils.estim.bounds_transform(10,inf),abs(utr(tr(11))-11)<1e-13
% [tr,utr]=utils.estim.bounds_transform(-inf,15),abs(utr(tr(11))-11)<1e-13

i0=isinf(lb) & isinf(ub);

i1=lb==0 & isinf(ub);

i2=isinf(lb) & ub==0;

i3=isfinite(lb) & isfinite(ub);

i4=~i1 & isfinite(lb) & isinf(ub);

i5=~i2 & isinf(lb) & isfinite(ub);

bad=~i0 & ~i1 & ~i2 & ~i3 & ~i4 & ~i5;

if any(bad)
    
    disp([lb(bad),ub(bad)])
    
    error('the constraints above are not supported')
    
end

is_i1=any(i1);

is_i2=any(i2);

is_i3=any(i3);

is_i4=any(i4);

is_i5=any(i5);

if all(i0)
    % no constraints
    tr=@(x)x;
    
    utr=@(x)x;
    
else
    
    tr=@transform;
    
    utr=@untransform;
    
end

    function x=transform(x)
        
        x(i1)=log(x(i1));
        
        x(i2)=log(-x(i2));
        
        x(i3)=log((x(i3)-lb(i3))./(ub(i3)-x(i3)));
        
        x(i4)=log(x(i4)-lb(i4));
        
        x(i5)=log(-x(i5)+ub(i5));
        
    end

    function x=untransform(x)
        
        if is_i1, x(i1)=exp(x(i1)); end
        
        if is_i2, x(i2)=-exp(x(i2)); end
        
        if is_i3
            
            ex3=exp(-x(i3));
            
            x(i3)=(ub(i3)+lb(i3).*ex3)./(1+ex3);
            
        end
        
        if is_i4, x(i4)=exp(x(i4))+lb(i4); end
        
        if is_i5, x(i5)=-exp(x(i5))+ub(i5); end
        
    end

end