function [keep,expand,C,UC,B]=shrink_expand(n,k,strategy,debug)
% shrink_expand computes shrinking and expansion objects for the manipulation of symmetric tensors
%
% Syntax
% -------
% ::
%
%   [keep,expand,C,UC,B]=shrink_expand(n,k)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy,debug)
%
% Inputs
% -------
%
% - **n** [scalar] : number of variables in the tensor
% - **k** [scalar] : order of the tensor
% - **strategy** [1|2|{3}] : alternative computation strategies
%   - 1: uses bsxfun
%   - 2: splanar inspired
%   - 3: applies ismember
% - **debug** [true|{false}] : checks the results
%
% Outputs
% --------
%
% - **keep** [logical] : n^k x 1 vector true for the columns to be kept
% - **expand** [vector] : n^k x 1 vector replicating the compressed
%   columns to form the grand tensor
% - **C** [matrix] : sparse compression matrix of size n^k x g, where
%   g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%   (matrix version of **keep**)
% - **UC** [matrix] : sparse expansion matrix of size g x n^k, where
%   g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%   (matrix version of **expand**)
% - **B** [matrix] : g x k matrix of combinations without repetitions. Each
%   row in increasing order.
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
            A=(1:nbig)*C;
            max(abs(A(:)-find(keep)))
        end
        
        if nargout>3
            % uncompression matrix
            %----------------------
            UC=speye(nkept);
            UC=UC(:,expand);
            if debug
                A=(1:nkept)*UC;
                max(abs(A(:)-expand(:)))
            end
            if nargout>4
                B=A(keep,:);
            end
        end
    end
end

    function ismember_strategy()
        [expand] = myismember(test2,test2(keep,:)); %<---[~,expand] = ismember(test2,test2(keep,:),'rows');
        function [locb] = myismember(A,B)
            % unique A first
            [icA] = myunique(); %<---[~,~,icA] = unique(A,'rows','sorted');
            
            % Sort the unique elements of B and B, duplicate entries are adjacent
            [sort_B_B,tags_B_B] = sortrows([B;B]);
            
            % Find matching entries
            d = sort_B_B(1:end-1,:)==sort_B_B(2:end,:);     
            d = all(d,2);                                   
            ndx1 = tags_B_B(d);                          
            
            % Find locb by using given indices
            locb = builtin('_ismemberfirst',icA,ndx1);
            function [indC] = myunique
                numRows = size(A,1);
                [sortA,indSortA] = sortrows(A);
                % groupsSortA indicates the location of non-matching entries.
                groupsSortA = sortA(1:numRows-1,:) ~= sortA(2:numRows,:);
                groupsSortA = any(groupsSortA,2);
                groupsSortA = [true; groupsSortA];          
                groupsSortA = full(groupsSortA);
                indC = cumsum(groupsSortA); 
                indC(indSortA) = indC; 
            end
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
