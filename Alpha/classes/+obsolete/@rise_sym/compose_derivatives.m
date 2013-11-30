function d=compose_derivatives(x,d_args)
        
switch x.func
    case 'plus'
        d=d_args{1}+d_args{2};
    case 'minus'
        d=d_args{1}-d_args{2};
    case 'mtimes'
        t11 = d_args{1}*x.args{2};
        t12 = d_args{2}*x.args{1};
        d=t11+t12;
    case 'mrdivide'
        if ~isa(d_args{2},'rise_sym')
            d_args{2}=rise_sym(d_args{2});
        end
        if ~is_zero(d_args{2})
            t11 = d_args{1}*x.args{2};
            t12 = d_args{2}*x.args{1};
            t13 = t11-t12;
            t14 = x.args{2}^2;
            d=t13/t14;
        else
            d=d_args{1}/x.args{2};
        end
    case {'lt','gt','le','ge','eq','ne','or','and'}
        d=rise_sym(0);
    case 'if_then_else'
        d=if_then_else(x.args{1},d_args{2},d_args{3});
    case 'if_elseif'
        the_args=x.args;
        the_args(2:2:end)=d_args(2:2:end);
        d=if_elseif(the_args{:});
    case 'mpower'
        t11 = log(x.args{1});
        t12 = d_args{2}*t11;
        t13 = d_args{1}*x.args{2};
        d= t12*x+t13*x.args{1}^(x.args{2}-1);
%         t14 = t13/x.args{1}; % <-- creates problems when x.args{1}=0
%         t15 = t12+t14;
%         d= t15*x.args{1}^x.args{2};
    case 'max'
        t11 = x.args{1}>x.args{2};
        t12 = t11*d_args{1};
        t13 = 1-t11;
        t14 = t13*d_args{2};
        d= t14+t12;
    case 'min'
        t11 = x.args{2}>x.args{1};
        t12 = t11*d_args{1};
        t13 = 1-t11;
        t14 = t13*d_args{2};
        d= t14+t12;
    case 'uminus'
        d= -d_args{1};
    case 'uplus'
        d= d_args{1};
    case 'exp'
        d= d_args{1}*x;
    case 'log'
        d= d_args{1}/x.args{1};
    case 'log10'
        t11 = exp(1);
        t12 = log10(t11);
        t13 = d_args{1}/x.args{1};
        d= t12*t13;
    case 'cos'
        t11 = sin(x.args{1});
        t12 = -t11;
        d= d_args{1}*t12;
    case 'sin'
        t11 = cos(x.args{1});
        d= d_args{1}*t11;
    case 'tan'
        t11 = x^2;
        t12 = 1+t11;
        d= d_args{1}*t12;
    case 'acos'
        t11 = sin(x);
        t12 = d_args{1}/t11;
        d= -t12;
    case 'asin'
        t11 = cos(x);
        d= d_args{1}/t11;
    case 'atan'
        t11 = x.args{1}^2;
        t12 = 1+t11;
        d= d_args{1}/t12;
    case 'cosh'
        t11 = sinh(x.args{1});
        d= d_args{1}*t11;
    case 'sinh'
        t11 = cosh(x.args{1});
        d= d_args{1}*t11;
    case 'tanh'
        d= d_args{1}*(1-x^2);
    case 'acosh'
        t11 = sinh(x);
        d= d_args{1}/t11;
    case 'asinh'
        t11 = cosh(x);
        d= d_args{1}/t11;
    case 'atanh'
        t11 = x.args{1}^2;
        t12 = 1-t11;
        d= d_args{1}*t12;
    case 'sqrt'
        t11 = 2*x;
        d= d_args{1}/t11;
    case 'abs'
        t11 = sign(x.args{1});
        d= t11*d_args{1};
    case 'sign'
        d= rise_sym(0);
    case 'erf'
        % x^2
        t11 = mpower(x.args{1},2);
        % exp(x^2)
        t12 =  exp(t11);
        % sqrt(pi)
        t11 = sqrt(pi);
        % sqrt(pi)*exp(x^2)
        t13 = t11*t12;
        % 2/(sqrt(pi)*exp(x^2));
        t14 = 2/t13;
        % (2/(sqrt(pi)*exp(x^2)))*dx;
        d= t14*d_args{1};
    case 'normcdf'
        t14 = sqrt(2*pi);
        % x - mu
        t12 = x.args{1}-x.args{2};
        % y = (x-mu)/sigma
        y = t12/x.args{3};
        % (x-mu)^2/sigma^2
        t12 = y^2;
        % -(x-mu)^2/sigma^2
        t13 = -t12;
        % -((x-mu)^2/sigma^2)/2
        t12 = t13/2;
        % exp(-((x-mu)^2/sigma^2)/2)
        t13 = exp(t12);
        % derivative of x.args{1} standardized normal
        % t15 = (1/sqrt(2*pi))*exp(-y^2/2)
        t15 = t13/t14;
        % derivatives thru x
        t11 = d_args{1}/x.args{3};
        % derivatives thru mu
        t12 = d_args{2}/x.args{3};
        % intermediary sum
        t14 = t11-t12;
        % derivatives thru sigma
        t11 = y/x.args{3};
        t12 = t11*d_args{3};
        % intermediary sum
        t11 = t14-t12;
        % total derivative:
        % (d_args{1}/sigma - d_args{2}/sigma - d_args{3}*(x-mu)/sigma^2) * t15
        % where t15 is the derivative of x.args{1} standardized normal
        d= t11*t15;
    case 'normpdf'
        % (x - mu)
        t11 = x.args{1}-x.args{2};
        % (x - mu)/sigma
        t12 = t11/x.args{3};
        % d_args{3} * (x - mu)/sigma
        t11 = d_args{3}*t12;
        % d_args{2} - d_args{1}
        t13 = d_args{2}-d_args{1};
        % d_args{2} - d_args{1} + d_args{3} * (x - mu)/sigma
        t14 = t13+t11;
        % ((x - mu)/sigma) * (d_args{2} - d_args{1} + d_args{3} * (x - mu)/sigma)
        t11 = t12*t14;
        % ((x - mu)/sigma) * (d_args{2} - d_args{1} + d_args{3} * (x - mu)/sigma) - d_args{3}
        t12 = t11-d_args{3};
        % this / sigma
        t11 = this/x.args{3};
        % total derivative:
        % (this / sigma) * (((x - mu)/sigma) * (d_args{2} - d_args{1} + d_args{3} * (x - mu)/sigma) - d_args{3})
        d= t11*t12;
end
% ensure that we remain rise_sym
if ~isa(d,'rise_sym')
    d=rise_sym(d); %d=rise_sym(d);
end

end
