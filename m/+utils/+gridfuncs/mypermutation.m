function p=mypermutation(v,p)
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

% finds all permutations of a vector v
% example v=[1,2,3];
% res=cell2mat(utils.gridfuncs.mypermutation(1:3));

if nargin<2

    p={[]};

end

n=numel(v);

np=numel(p);

if n==1
    
    for ip=1:np
    
        p{ip}=[p{ip},v];

    end
    
else
    
    iter=0;
    
    for ii=1:n
        
        vi=v;
        
        vi(ii)=[];
        
        ppi=utils.gridfuncs.mypermutation(vi,p);
        
        if ii==1
            
            nppi=numel(ppi);
            
            pp_nbr=nppi*n;
        
            pp=cell(pp_nbr,1);
        
        end
        
        for jj=1:nppi
            
            iter=iter+1;
        
            pp{iter}=[v(ii),ppi{jj}];
    
        end
        
    end
    
    p=pp;

end

end