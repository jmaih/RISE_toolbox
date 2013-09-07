clear classes
clc
%% choose level
level=1;
switch level
    case 1
        funcs={'sin(x)*exp(cos(log(y)))*z+k^3*(l-x)^2';'atan(x)*exp(cos(log(y)))*z'};
        varList={'x','y','z','k','l'};
        wrt={'z','y','x','k','l'};
    case 2
        load difftest
        funcs=symbolic_original;
        varList=symb_list;
        wrt=with_respect_to;
    case 3
        load difftest2
        funcs=symbolic_batch;
        varList=symb_list;
        wrt=with_respect_to;
end
eqtn_nbr=numel(funcs);
wrt_nbr=numel(wrt);

profile off
profile on
args=strcat(varList,',');
args=cell2mat(args);

tic
[mytree,wrt_index,wrt_,mycall]=sadiff.setup(funcs,varList,wrt);
t0=toc;
fprintf(1,'%s %0.10g\n','tree construction ',t0);

tic
this=diff(mytree,wrt_);
t1=toc;
fprintf(1,'%s %0.10g\n','reverse differentiation of the tree ',t1);

tic
[c,mycall]=print(this);
t2=toc;
JJ1=strcat(mycall.fid(:,1),'=',mycall.fid(:,2),';');
for idef=1:numel(c)
    JJ1=strrep(JJ1,c{idef},['Jac(',sprintf('%0.10g',idef),',indx)']);
end
JJ1=[['Jac=zeros(',int2str(eqtn_nbr),',',int2str(wrt_nbr),');']
    ['bigi_=speye(',int2str(wrt_nbr),');']
    JJ1];
fprintf(1,'%s %0.10g\n','printing of reverse diff ',t2);
fprintf(1,'%s %0.10g\n','reverse diff + printing of tree ',t1+t2);

tic
[xhandle,dxhandle,~,~,mycall] = differentiate(mytree,wrt_,wrt_index);
t3=toc;
JJ2=strcat(mycall.fid(:,1),'=',mycall.fid(:,2),';');
for idef=1:numel(dxhandle)
    JJ2=strrep(JJ2,dxhandle{idef},['Jac(',sprintf('%0.10g',idef),',indx_',int2str(idef),'_)']);
end
JJ2=[['Jac=zeros(',int2str(eqtn_nbr),',',int2str(wrt_nbr),');']
    ['bigi_=speye(',int2str(wrt_nbr),');']
    JJ2];
fprintf(1,'%s %0.10g\n','forward differentiation of tree and printing ',t3);

profile off
profile viewer

%% expansion of derivatives in case there are many wrt and if we want to go to higher order...
% this function cannot expand many objects simultaneously. it stops when it
% finds a vector of objects. this is why we loop
tic
tmp=sadiff.empty(0);
for ii=1:numel(this)
    tmp(ii,:)=expand(this(ii));
end
toc

%% print expanded version of the jacobian (without the intermediate steps)

JJ0=char(this,[],[],true);

%% further examples

% raw
[jaco_tree3,code3,code_expanded3]=sadiff.jacobian(funcs,varList,wrt);

% tree and wrt done
[jaco_tree4,code4,code_expanded4]=sadiff.jacobian(mytree,varList,wrt_);

% hybrid (tree not done and wrt done)
[jaco_tree5,code5,code_expanded5]=sadiff.jacobian(funcs,varList,wrt_);

% other way around
[jaco_tree6,code6,code_expanded6]=sadiff.jacobian(mytree,varList,wrt);


