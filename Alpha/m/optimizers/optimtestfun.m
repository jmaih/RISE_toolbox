function y=optimtestfun(x,choice,shift)
error(nargchk(1,3,nargin))
data_flag=false;
if ischar(x)
    n=1;
else
    n=numel(x);
end
if nargin<3
    shift=[];
    if nargin<2
        choice=x;
     data_flag=true;
   end
end

if ~data_flag
    if isempty(shift)
        shift=0;
    end
    x=x-shift;
end
switch lower(choice)
    case {1,'sphere'}
        if data_flag
            y={'Sphere',100*[-1,1]};
        else
            y=sum(x.^2);
        end
    case {2,'f2'}
        if data_flag
            y={'f2',10*[-1,1]};
        else
            ax=abs(x);
            y=(sum(ax)+prod(ax));
        end
    case {3,'ridge'}
        if data_flag
            y={'Ridge',100*[-1,1]};
        else
            y=0;
            for ii=1:n
                y=y+sum(x(1:ii))^2;
            end
        end
    case {4,'f4'}
        if data_flag
            y={'f4',100*[-1,1]};
        else
            y=max(abs(x));
        end
    case {5,'rosenbrock'}
        if data_flag
            y={'Rosenbrock',30*[-1,1]};
        else
            x1=x(1:end-1)-1;
            xx=x(2:end)-x(1:end-1).^2;
            y=sum(100*xx.^2+x1.^2);
        end
    case {6,'f6'}
        if data_flag
            y={'f6',100*[-1,1]};
        else
            y=sum(floor(x+.5).^2);
        end
    case {7,'f7'}
        if data_flag
            y={'f7',128*[-1,1]};
        else
            y=sum((1:n)'.*x(:).^4);
        end
    case {8,'schwefel'}
        if data_flag
            y={'Schwefel',500*[-1,1]};
        else
            y=-sum(x.*sin(sqrt(abs(x))));
        end
    case {9,'rastrigin'}
        if data_flag
            y={'Rastrigin',5.12*[-1,1]};
        else
            y=sum(x.^2-10*cos(2*pi*x)+10);
        end
    case {10,'ackley'}
        if data_flag
            y={'Ackley',32*[-1,1]};
        else
            y=20-20*exp(-.2*sqrt(sum(x.^2)/n))+exp(1)-exp(sum(cos(2*pi*x))/n);
        end
    case {11,'griewank'}
        if data_flag
            y={'Griewank',600*[-1,1]};
        else
            y=sum(x.^2)/4000-prod(cos(x(:)./sqrt((1:n)')))+1;
        end
    case {12,'levymontalvo1'}
        if data_flag
            y={'LevyMontalvo1',10*[-1,1]};
        else
            yi=1+.25*(x+1);
            py=pi*yi;
            y=pi/n*(...
                10*sin(py(1)).^2+...
                sum(...
                (yi(1:end-1)-1).^2.*(1+10*sin(py(2:end)).^2)...
                )+...
                (yi(end)-1)^2);
        end
    case {13,'levymontalvo2'}
        if data_flag
            y={'LevyMontalvo2',5*[-1,1]};
        else
            px=pi*x;
            y=0.1*(...
                10*sin(3*px(1)).^2+...
                sum(...
                (x(1:end-1)-1).^2.*(1+sin(3*px(2:end)).^2)...
                )+...
                (x(end)-1)*(1+sin(2*px(1))^2));
        end
    case {14,'f14'}
        if data_flag
            y={'f14',100*[-1,1]};
        else
            y=sum(x.^2)+(sum(.5*(1:n)'.*x(:)))^2+...
                (sum(.5*(1:n)'.*x(:)))^4;
        end
    case {15,'cosinemixture'}
        if data_flag
            y={'CosineMixture',1*[-1,1]};
        else
            y=sum(x.^2)-sum(cos(5*pi*x));
        end
    case {16,'epistaticmichalewicz'}
        if data_flag
            y={'EpistaticMichalewicz',1*[0,pi]};
        else
			theta=pi/6;
			m=10;
			ii=(1:n)';
			odd=rem(ii(1:end-1),2)~=0;
			even=~odd;
			yio=x(1:end-1)*cos(theta)-x(2:end)*sin(theta);
			yie=x(1:end-1)*cos(theta)+x(2:end)*sin(theta);
			yi=yio;
			yi(even)=yie(even);
			yi=[yi;x(end)];
            y=-sum(sin(yi).*(sin(ii.*yi.^2/pi)).^(2*m));
        end
    case {17,'exponential'}
        if data_flag
            y={'Exponential',1*[-1,1]};
        else
            y=-exp(-.5*sum(x.^2));
        end
    case {18,'salomon'}
        if data_flag
            y={'Salomon',100*[-1,1]};
        else
			nx=norm(x);
            y=1-cos(2*pi*nx)+.1*nx;
        end
    case {19,'shubert'}
        if data_flag
            y={'Shubert',10*[-1,1]};
        else
			jj=1:5;
			y=1;
			for ii=1:n
				y=y*(sum(jj.*cos((jj+1)*x(ii)+jj)));
			end
        end
    otherwise
end

% function u=u_func(x,a,k,m)
% u=0*x;
% h=x>a;
% l=x<-a;
% u(h)=k*(x(h)-a).^m;
% u(l)=k*(-x(l)-a).^m;