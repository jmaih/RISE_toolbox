function d=diff(x,wrt,alien_list)
% INTERNAL FUNCTION: overloads diff for splanar
%

if nargin<3

    alien_list=[];

end


[nrows,ncols]=size(x);

if nrows>1||ncols>1

    d=cell(nrows,ncols);

    for irow=1:nrows

        for icol=1:ncols

            d{irow,icol}=diff(x(irow,icol),wrt,alien_list);

        end

    end

    return

end

nargs=numel(x.args);

if isempty(x.incidence)
    % numbers/vectors and variables which are not part of differentiation
    % automatically receive 0 as derivative
    d=splanar.prototypize();

elseif nargs==0
    % variables that are part of differentiation
    d=splanar.prototypize(double(x.incidence(wrt)));
    % why do we need to double it?
    % we need to carry around vectors in order to be able to do
    % higher-order derivatives

elseif ~isempty(alien_list) && ismember(x.func,alien_list)

    d=x;

    is_differentiated=false;

    for ii=1:nargs

        d_arg=d.args{ii};

        is_differentiated=ischar(d_arg) && strcmp(d_arg,'''diff''');

        if is_differentiated

            break

        end

    end

    if is_differentiated
        % increase the order of differentiation
        d.args{end}=d.args{end}+1;

    else
        % first-order derivative
        d.args=[d.args,{'''diff''',1}];

    end

else
    % functions of variables
    if_elseif_flag=strcmp(x.func,'if_elseif');

    if_then_else_flag=strcmp(x.func,'if_then_else');

    d_args=x.args;

    for iarg=1:nargs

        if (if_elseif_flag && rem(iarg,2))||(iarg==1 && if_then_else_flag)

            continue

        end

        d_args{iarg}=diff(x.args{iarg},wrt,alien_list);

    end
    % compose derivatives
    %--------------------
    d=compose_derivative();%x,d_args

end

    function d=compose_derivative()%x,d_args

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

                if ~is_zero(d_args{2})

                    t11 = d_args{1}*x.args{2};

                    t12 = d_args{2}*x.args{1};

                    t13 = t11-t12;

                    t14 = x.args{2}^2;

                    d=t13/t14;

                else

                    d=d_args{1}/x.args{2};

                end

            case {'lt','gt','le','ge','eq','ne','or','and','sign',...
                    'steady_state'}

                % derivatives are zero
                d=splanar.prototypize();

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

                d= d_args{1}*(1+x^2);

            case 'cot'

                d= -d_args{1}*(1-x^2); % <--- d=-d_args{1}/sin(x.args{1})^2;

            case 'acos'

                d= -d_args{1}/sqrt(1-x^2);

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

            case 'norminv'

                d=(d_args{1}+d_args{2}+d_args{3})/normpdf(...
                    norminv(x.args{1},x.args{2},x.args{3}));

            case 'normcdf'

                d=(d_args{1}+d_args{2}+d_args{3})*normpdf(x.args{1},x.args{2},x.args{3});

            case 'normpdf'

                y=(x.args{1}-x.args{2})/x.args{3};

                ss= d_args{3}/x.args{3};

                d=((ss*y-(d_args{1}-d_args{2})/x.args{3})*y-ss)*x;

            case 'betainv'

                d=(d_args{1}+d_args{2}+d_args{3})/betapdf(...
                    betainv(x.args{1},x.args{2},x.args{3}));

            case 'betacdf'

                d=(d_args{1}+d_args{2}+d_args{3})*betapdf(x.args{1},x.args{2},x.args{3});

            case 'betapdf'

                numerator=(x.args{2}+x.args{3}-2)*x.args{1}+1-x.args{2};

                denominator= (x.args{1}-1)*x.args{1};

                d=numerator/denominator*betapdf(x.args{1},x.args{2},x.args{3});

        end

    end

end