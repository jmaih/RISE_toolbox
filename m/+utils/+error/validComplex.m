function flag=validComplex(x)
% INTERNAL FUNCTION: check validity of values of different data types
%

% This function is the same as utils.error.valid except that it allows
% complex numbers

check_real=false;

flag=utils.error.valid(x,check_real);

end