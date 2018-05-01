function lhs=korder_matrix_vector(x,T_order_1,G,G0p)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


[r1,c1,p1]=size(T_order_1);
[rg,cg,pg,bg]=size(G);
[rg0p,cg0p,pg0p]=size(G0p);

h=bg;
if (h~=pg)||(h~=p1)||(h~=pg0p)||...
        (cg0p~=cg)||(rg0p~=rg)||(r1~=rg)||...
        (r1~=c1)
    error('wrong matrices formats')
end

ru=cg;
nx=numel(x);
cu=nx/(ru*h);
k=cu/r1;

U=reshape(x,[ru,cu,h]); % rename U to x in order to save on memory...

lhs=zeros(ru,cu*h);
for st=1:h
    tmp=0;
    for stp=1:h
       tmp=tmp+G(:,:,st,stp)*U(:,:,stp);
    end
    tmp=utils.kronecker.A_times_k_kron_B(tmp,T_order_1(:,:,st),k);
    tmp=tmp+G0p(:,:,st)*U(:,:,st);
    lhs(:,(st-1)*cu+1:st*cu)=tmp;
end

lhs=lhs(:);