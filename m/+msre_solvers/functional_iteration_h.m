function [T1,T0_T1]=functional_iteration_h(T0,dbf_plus,d0,dpb_minus,bf_cols,pb_cols,~)
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


n=size(d0{1},1);

h=size(d0,2);

if nargin<6
    
    pb_cols=[];
    
    if nargin<5
        
        bf_cols=[];
        
    end
    
end

if isempty(bf_cols)
    
    bf_cols=1:n;
    
end

if isempty(pb_cols)
    
    pb_cols=1:n;
    
end

npb=numel(pb_cols);

if size(dbf_plus{1},2)~=numel(bf_cols)
    
    error('number of columns of dbf_plus inconsistent with the number of bf variables')
    
end

if size(dpb_minus{1},2)~=npb
    
    error('number of columns of dpb_minus inconsistent with the number of bp variables')
    
end

T0=reshape(T0,[n,npb,h]);

T1=T0;

for r0=1:h
    
    U=d0{r0};
    
    for r1=1:h
        
        U(:,pb_cols)=U(:,pb_cols)+dbf_plus{r0,r1}*T0(bf_cols,:,r1);
        
    end
    
%     T1(:,:,r0)=-pinv(full(U))*dpb_minus{r0}; 
    T1(:,:,r0)=-U\dpb_minus{r0};
    
end

% update T
%---------
T1=reshape(T1,[n,npb*h]);

T0_T1=reshape(T0,[n,npb*h])-T1;

end