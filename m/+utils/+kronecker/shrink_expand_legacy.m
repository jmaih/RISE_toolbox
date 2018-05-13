function [keep,expand,C,UC,B]=shrink_expand_legacy(n,k,strategy,debug)
% shrink_expand computes shrinking and expansion objects for the manipulation of symmetric tensors
%
% ::
%
%
%   [keep,expand,C,UC,B]=shrink_expand(n,k)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy,debug)
%
% Args:
%
%    - **n** [scalar] : number of variables in the tensor
%    - **k** [scalar] : order of the tensor
%    - **strategy** [1|2|{3}] : alternative computation strategies
%      - 1: uses bsxfun
%      - 2: splanar inspired
%      - 3: applies ismember
%    - **debug** [true|{false}] : checks the results
%
% Returns:
%    :
%
%    - **keep** [logical] : n^k x 1 vector true for the columns to be kept
%    - **expand** [vector] : n^k x 1 vector replicating the compressed
%      columns to form the grand tensor
%    - **C** [matrix] : sparse compression matrix of size n^k x g, where
%      g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%      (matrix version of **keep**)
%    - **UC** [matrix] : sparse expansion matrix of size g x n^k, where
%      g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%      (matrix version of **expand**)
%    - **B** [matrix] : g x k matrix of combinations without repetitions. Each
%      row in increasing order.
%
% Note:
%
%    useful for shrinking tensors of the form fvvv...v as used in higher-order
%    differentiation
%
% Example:
%
%    See also:

if nargin<4
    
    debug=false;
    
    if nargin<3
        
        strategy=3;
        
    end
    
end

% differentiation order all derivatives
%--------------------------------------
A=utils.gridfuncs.mygrid(n*ones(1,k));

nbig=size(A,1);

% find the non-decreasing indices
%--------------------------------
% flip left right
test2=A(:,end:-1:1);
% stamp the decreasing rows
drow=test2(:,1:end-1)-test2(:,2:end);

% definition change to accord with derivatives
%----------------------------------------------
keep=all(drow<=0,2);%<--keep=~any(drow<0,2);

if nargout>1
    
    nkept=sum(keep);
    
    % now map all rows into the kept ones and vice versa
    %----------------------------------------------------
    test2=sort(test2,2);
    
    % separate kept and unkept
    kept=test2(keep,:);
    
    % stamps
    
    switch strategy
        
        case 1
            
            expand=bsxfun_strategy();
            
        case 2
            
            expand=splanar_strategy();
            
        case 3
            
            [~,expand]=utils.gridfuncs.ismember(test2,kept);
            
        otherwise
            
            error('strategy not implemented')
            
    end
    
    if nargout>2
        % compression matrix
        %-------------------
        C=speye(nbig);
        
        C=C(:,keep);
        
        if debug
            
            A_=(1:nbig)*C;
            
            max(abs(A_(:)-find(keep)))
            
        end
        
        if nargout>3
            % uncompression matrix
            %----------------------
            UC=speye(nkept);
            
            UC=UC(:,expand);
            
            if debug
                
                A_=(1:nkept)*UC;
                
                max(abs(A_(:)-expand(:)))
                
            end
            
            if nargout>4
                
                B=A(keep,:);
                
            end
            
        end
        
    end
    
end

    function expand=splanar_strategy()
        
        expand=zeros(1,nbig);
        
        proto_permutation=cell2mat(utils.gridfuncs.mypermutation(1:k));
        
        for ikept=1:nkept
            
            this=kept(ikept,:);
            
            this=this(proto_permutation);
            
            % locate these permutations
            ypred=utils.gridfuncs.locate_permutation(this,n,false);
            
            ypred=unique(ypred);
            
            expand(ypred)=ikept;
            
        end
        
    end

    function expand=bsxfun_strategy()
        
        expand=zeros(1,nbig);
        
        expand(keep)=1:nkept;
        
        Index=1:nbig;
        
        Index(keep)=[];
        
        if ~isempty(Index)
            
            for ikept=1:nkept
                
                target=kept(ikept,:);
                
                if all(target==target(1))
                    
                    continue
                    
                end
                
                bingo=sum(abs(bsxfun(@minus,test2(Index,:),target)),2)==0;
                
                expand(Index(bingo))=ikept;
                
                Index(bingo)=[];
                
            end
            
        end
        
    end

end
% not implemented solution : create a container with unique identifiers
