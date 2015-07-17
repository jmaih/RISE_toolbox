function val=polynomial_integration(p,a,b)
% polynomial is of the form a0+a1*x+...+ar*x^r
% the integral is then a0*x+a1/2*x^2+...+ar/(r+1)*x^(r+1)

order=numel(p)-1;

pbar=p(:)./(1:order+1)';

b_a=b-a;
val=0;
% exploit horner
for oo=order+1:-1:1
    val=(val+pbar(oo))*b_a;
end

