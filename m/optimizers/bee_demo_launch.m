clear
clc
dbclear if error
close all
s = RandStream.create('mt19937ar','seed',1974);
RandStream.setDefaultStream(s);

funcs={'schwefel','rastrigin','ackley','griewank'};

n=2; % number of parameters
opt=struct();
opt.colony_size=10;
% opt.B=10;
opt.max_fcount=inf;
opt.max_iter=100;
opt.display=1;
opt.fbest=0;

func=funcs{2};
datta=optimtestfun(func);
lb=datta{2}(1)*ones(n,1);
ub=datta{2}(2)*ones(n,1);
shift=(lb+(ub-lb).*rand(n,1));
lb=lb+shift;
ub=ub+shift;
opt.xbest=shift;

x0=lb+(ub-lb).*rand(n,1);
opt.griddata=cell(1,n+1);
for ii=1:n
    opt.griddata{ii}=linspace(lb(ii),ub(ii),round(40/n));
    opt.griddata{ii}=[opt.griddata{ii},shift(ii)];
    opt.griddata{ii}=sort(opt.griddata{ii},2,'ascend');
end
if n==2
    [opt.griddata{1},opt.griddata{2}]=meshgrid(opt.griddata{1},opt.griddata{2});
end
opt.griddata{end}=opt.griddata{end-1};
for ii=1:numel(opt.griddata{end-1});
    xi=[];
    for jj=1:n
        xi=[xi;opt.griddata{jj}(ii)];
    end
    opt.griddata{end}(ii)=optimtestfun(xi,func,shift);
end
% dhms(@optimtestfun,x0,lb,ub,opt,func);

obj=bee_demo(@optimtestfun,x0,[],lb,ub,opt,func,shift);
