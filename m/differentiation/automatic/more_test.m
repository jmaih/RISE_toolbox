clc
x=rand(10,1);
dx=eye(numel(x));
occur=[1,3,5,6];
testfun=@(x)x(1)*log(x(3))-x(5)*exp(x(6));
aa=automatic(x);
testfun2=@(x)x(occur(1))*log(x(occur(2)))-x(occur(3))*exp(x(occur(4)));
testfun3=@(x)x(1)*log(x(2))-x(3)*exp(x(4));

t0=[0,0,0];
neq=300;
for ii=1:neq
tic,a_result=testfun(aa);t0(1)=toc;
tic,bb=automatic(x,dx(:,occur));b_result=testfun2(bb);t0(2)=toc;
tic,cc=automatic(x(occur),dx(occur,occur)); c_result=testfun3(cc);t0(3)=toc;
end
disp(a_result)
disp(b_result)
disp(c_result)
disp(t0/min(t0))