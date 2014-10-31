function [keep,expand,C,UC]=shrink_expand(n,k,strategy,debug)
% shrink_expand computes shrinking and expansion objects for the manipulation of symmetric tensors
%
% Syntax
% -------
% ::
%
%   [keep,expand,C,UC]=shrink_expand(n,k)
%   [keep,expand,C,UC]=shrink_expand(n,k,strategy)
%   [keep,expand,C,UC]=shrink_expand(n,k,strategy,debug)
%
% Inputs
% -------
%
% - **n** [scalar] : number of variables in the tensor
% - **k** [scalar] : order of the tensor
% - **strategy** [1|2|3|4|{5}] : alternative computation strategies
%   - 1: bsxfun
%   - 2: splanar inspired
%   - 3: cumul
%   - 4: cumul
%   - 5: ismember
% - **debug** [true|{false}] : checks the results
%
% Outputs
% --------
%
% - **keep** [logical] : n^k x 1 vector true for the columns to be kept
% - **expand** [vector] : n^k x 1 vector replicating the compressed
%   columns to form the grand tensor
% - **C** [matrix] : sparse compression matrix of size n^k x g, where g is
%   the number of unique elements in the tensor (matrix version of
%   **keep**)
% - **UC** [matrix] : sparse expansion matrix of size g x n^k g, where g is
%   the number of unique elements in the tensor (matrix version of
%   **expand**)
%
% More About
% ------------
%
% useful for shrinking tensors of the form fvvv...v as used in higher-order
% differentiation
%
% Examples
% ---------
%
% See also:


if nargin<4
    debug=false;
    if nargin<3
        strategy=5;
    end
end

% differentiation order all derivatives
%--------------------------------------
test=utils.gridfuncs.mygrid(n*ones(1,k));
nbig=size(test,1);

% find the non-decreasing indices
%--------------------------------
% flip left right
test2=test(:,end:-1:1);
% stamp the decreasing rows
drow=test2(:,1:end-1)-test2(:,2:end);
keep=~any(drow<0,2);
if nargout>1
    nkept=sum(keep);
    
    % now map all rows into the kept ones and vice versa
    %----------------------------------------------------
    test2=sort(test2,2);
    
    % stamps
    expand=nan(1,nbig);
    
    % separate kept and unkept
    kept=test2(keep,:);
    
    switch strategy
        case 1
            bsxfun_strategy()
        case 2
            splanar_strategy()
        case 3
            cumul_strategy()
        case 4
            cumul_strategy2()
        case 5
            ismember_strategy()
        otherwise
            error('strategy not implemented')
    end
    
    if nargout>2
        % compression matrix
        %-------------------
        C=speye(nbig);
        C=C(:,keep);
        if debug
            test=(1:nbig)*C;
            max(abs(test(:)-find(keep)))
        end
        
        if nargout>3
            % uncompression matrix
            %----------------------
            UC=speye(nkept);
            UC=UC(:,expand);
            if debug
                test=(1:nkept)*UC;
                max(abs(test(:)-expand(:)))
            end
        end
    end
end

    function ismember_strategy()
        A=sum(cumsum(test2,2),2);
        B=A(keep);
        [~,expand] = ismember(A,B);
    end

    function cumul_strategy2()
        tmp=sum(cumsum(test2,2),2);
        tmp_keep=tmp(keep);
        Index=1:nbig;
        Index(keep)=[];
        expand(keep)=1:nkept;
        for ikept=1:nkept
            match=tmp(Index)==tmp_keep(ikept);
            expand(Index(match))=ikept;
            Index(match)=[];
        end
    end

    function cumul_strategy()
        tmp=sum(cumsum(test2,2),2);
        tmp_keep=tmp(keep);
        for ikept=1:nkept
            expand(tmp==tmp_keep(ikept))=ikept;
        end
    end

    function splanar_strategy()
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

    function bsxfun_strategy()
        expand(keep)=1:nkept;
        Index=1:nbig;
        Index(keep)=[];
        for ikept=1:nkept
            bingo=sum(abs(bsxfun(@minus,test2(Index,:),kept(ikept,:))),2)==0;
            expand(Index(bingo))=ikept;
            Index(bingo)=[];
        end
    end
end
% not implemented solution : create a container with unique identifiers
