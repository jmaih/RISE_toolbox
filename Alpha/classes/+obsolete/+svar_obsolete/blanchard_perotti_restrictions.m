function [R,lb,ub]=blanchard_perotti_restrictions(param)

if nargin<1
    R={'APITR','GOV','RGDP'};
else
    a_y=2.08;
    a_g=param(1); % -inf inf
    sig_t=param(2); % 0 inf
    b_y=param(3); % -inf inf
    b_t=param(4); % -inf inf
    sig_g=param(5); % 0 inf
    c_t=param(6); % -inf inf
    c_g=param(7); % -inf inf
    sig_y=param(8); % 0 inf
    A=[1 0 -a_y
        0 1 -b_y
        -c_t -c_g 1];
    
    B=[sig_t sig_g*a_g 0
        b_t*sig_t sig_g 0
        0 0 sig_y];
    R=A\B;
end

if nargout>1
    ub=inf*ones(8,1);
    lb=-ub;
    lb([2,5,8])=sqrt(eps);
end


