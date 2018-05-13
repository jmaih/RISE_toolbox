function derivs=differentiate(eqtns,nwrt,order,verbose)
% differentiate - differentiates vectors of splanar objects
%
% ::
%
%
%   derivs=splanar.differentiate(eqtns,nwrt)
%   derivs=splanar.differentiate(eqtns,nwrt,order)
%   derivs=splanar.differentiate(eqtns,nwrt,order,verbose)
%
% Args:
%
%    - **eqtns** [splanar|cell array]: vector or cell array of splanar objects
%
%    - **nwrt** [integer]: number of variables for which differentiation is
%      taken
%
%    - **order** [integer|{1}]: order of differentiation
%
%    - **verbose** [true|{false}]: displays information about the process e.g.
%      the amount of time it takes to differentiate each order
%
% Returns:
%    :
%
%    - **derivs** [structure]: each element in the structure contains:
%      - **size** [2-colum vector]: number of rows and number of columns of
%        the compacted derivatives i.e. unique columns
%      - **derivatives** [vector of splanar]: derivatives
%      - **nwrt** [integer]: number of variables for which differentiation is
%        taken
%      - **order** [integer]: order of differentiation of the current
%        structure
%      - **nnz_derivs** [integer]: number of non-zero derivatives
%      - **partitions** [vector]: vector that help reconstruct the
%        expanded/grand derivative matrix.
%      - **map** []: empty field that will be used in a later stage
%
% Note:
%
% Example:
%
%    See also:

if nargin<4

    verbose=false;

end

if iscell(eqtns)

    eqtns=[eqtns{:}];

end

eqtns=eqtns(:);

neqtns0=size(eqtns,1);

% identify where the equations are coming from since everything will be
% vectorized later on
for ieq=1:neqtns0

    eqtns(ieq).location={ieq,[]};

end

neqtns=neqtns0;

derivs=struct('size',{},'derivatives',{},...
    'nwrt',{},'order',{},'nnz_derivs',{},'partitions',{},...
    'map',{},'vectorizer',{});

[keep_all,expand_all]=utils.kronecker.shrink_expand(nwrt,order);

for oo=1:order
    
    if verbose
    
        tic
    
    end
    
    % compute all combinations and hash them up
    %------------------------------------------
    [ncols,B,expand]=store_combinations();
    
    % differentiate
    %--------------
    d=cell(nwrt*neqtns,1);
    
    iter=0;
    
    proto_offset=0;
    
    nnz_derivs=0;
    
    big_index=cell(1,neqtns);
    
    for ieqtn=1:neqtns % number of rows
        
        oldwrt=eqtns(ieqtn).lineage;
        
        if ~isempty(oldwrt)
        
            oldwrt=oldwrt{end};
        
        end
        
        ncols_=max(1,numel(oldwrt));
        
        index_=cell(1,ncols_);
        
        for icol=1:ncols_ % number of derivatives taken in earlier round
            
            obj=intercept_column(eqtns(ieqtn),icol);
            
            if ~isempty(eqtns(ieqtn).lineage) && ~isempty(eqtns(ieqtn).lineage{end})
            
                obj.lineage=[eqtns(ieqtn).lineage(1:end-1),{eqtns(ieqtn).lineage{end}(icol)}];
            
            end
            
            candidates=find(obj.incidence);

            % ensure the location remains the same as in the parent
            % object. otherwise the equation number will be lost
            %-------------------------------------------------------
            obj.location=eqtns(ieqtn).location;
            
            wrt=candidates;
            
            if ~isempty(wrt)
                
                if ~isempty(oldwrt) % <--- oo>1
                    
                    wrt(wrt<oldwrt(icol))=[];
                
                    lineage__=[obj.lineage(1:end-1),{oldwrt(icol),wrt}];
                
                else
                    
                    lineage__={wrt};
                
                end
                
                d0=diff(obj,wrt);
                
                d0.lineage=lineage__;
                
                nderivs=numel(d0.lineage{end});
                % store the derivative if not zero
                %---------------------------------
                if nderivs
                    % locate each derivative in the reduced matrix
                    %---------------------------------------------
                    former=cell2mat(d0.lineage(1:end-1));
                    
                    if isempty(former)
                    
                        current=lineage__{end}(:);
                    
                    else
                        
                        current=[former(ones(nderivs,1),:),lineage__{end}(:)];
                    
                    end
                    
                    deriv_positions=proto_offset+(1:nderivs);
                    
                    index_{icol}=set_positions(current);
                    
                    d0.location={obj.location{1},deriv_positions};
                    
                    proto_offset=deriv_positions(end);
                    
                    iter=iter+1;
                    
                    d{iter}=d0;
                
                    nnz_derivs=nnz_derivs+nderivs;
            
                end
                
            end
            
        end
        
        big_index{ieqtn}=cell2mat(index_);
        
    end
    
    % make sure we have a row vector here!!!
    %----------------------------------------    
    big_index=cell2mat(big_index);
    
    d=[d{1:iter}];
    
    d=d(:);
    
    % swap locations
    %---------------
    neqtns=numel(d);
    
    for ideriv=1:neqtns
    
        d(ideriv).location{2}=big_index(d(ideriv).location{2});
    
    end
    
    % spit out the time it took to compute
    %-------------------------------------
    if verbose
    
        fprintf(1,'differentiation at order %0.0f done in %0.4f seconds\n',oo,toc);

    end
    
    derivs(oo)=struct('size',{[neqtns0,ncols]},'derivatives',d,...
        'nwrt',nwrt,'order',oo,'nnz_derivs',nnz_derivs,'partitions',expand,...
        'map',[],'vectorizer',[]);
    
    % next round
    %-----------
    eqtns=d;

end

    function p=set_positions(proto)
        
        n=size(proto,1);
        
        p=zeros(1,n);
        
        batch=num2cell(proto);
        
        for ii=1:n
            
            p(ii)=B(batch{ii,:});
            
        end
        
    end

    function [nkept,B,expand]=store_combinations()
        
        keep=keep_all{oo};
        
        expand=expand_all{oo};
        
        if oo==1
            
            B=expand(:);
            
        else
            
            tmp=repmat({nwrt},1,oo);
            
            B=reshape(expand,[tmp{:}]);
            
        end
    
        nkept=sum(keep); % size(B,1)

    end

end

% % we keep it outside the space of the main function otherwise memoization
% % will store useless and heavy too many variables in its workspace.
% function expand=inflator(nwrt,oo)
% [~,expand]=utils.kronecker.shrink_expand(nwrt,oo);
% end
