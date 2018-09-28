clc
nmat=4;
imax=10;
ms=randi(imax,nmat,2);
Ai=cell(1,nmat);
for ii=1:nmat
    Ai{ii}=rand(ms(ii,1),ms(ii,2));
end
newOrder=randperm(nmat);
k=@utils.kronecker.kronall;
A=k(Ai{:});
tic,B=utils.kronecker.shuffle_tensor1(A,ms,newOrder);toc
max(max(abs(B-k(Ai{newOrder}))))

sp=@utils.kronecker.sum_permutations;
opts=struct('skip_first',true);

opts.algo='shuf1';
tic,P1=sp(A,ms,opts,newOrder);toc
max(max(abs(B-P1)))

opts.algo='shuf2';
tic,P2=sp(A,ms,opts,newOrder);toc
max(max(abs(B-P2)))

opts.algo='old';
tic,P3=sp(A,ms,opts,newOrder);toc
max(max(abs(B-P3)))

opts.algo='perm';
tic,P4=sp(A,ms,opts,newOrder);toc
max(max(abs(B-P4)))


