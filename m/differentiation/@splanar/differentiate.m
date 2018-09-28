function derivs=differentiate(eqtns,nwrt,order,alien_list,verbose)
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
%    - **alien_list** [empty|cellstr|char]: list of alien functions
%      (functions that RISE does not recognize and that are to be
%      differentiated by the user himself).
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

% container by order and by equation and not container for all
% equations in all orders...
%
% We simultaneously count the number of nonzero derivatives
%
% we cannot pair the derivatives in the expanded matrix...
% since we are computing vectors belonging to different
% locations. Unless at evaluation, we create the vector and
% then dispatch each element separately, which can also be
% done.
%
% We can offer both options at the time of writing
% derivatives...

if nargin<5
    
    verbose=[];
    
    if nargin<4
        
        alien_list=[];
        
        if nargin<3
            
            order=[];
            
        end
        
    end
    
end

if isempty(verbose),verbose=false; end

if isempty(order),order=1; end

if iscell(eqtns)
    
    eqtns=[eqtns{:}];
    
end

eqtns=eqtns(:);

neqtns=size(eqtns,1);

% identify where the equations are coming from since everything will be
% vectorized later on
for ieq=1:neqtns
    
    eqtns(ieq).location={ieq,[]};
    
end

derivs=struct('size',{},'derivatives',{},...
    'nwrt',{},'order',{},'nnz_derivs',{},...  
    'timer',{});

perms=cell(1,order);

for ieqtn=1:neqtns
    
    myeqtn=eqtns(ieqtn);
    
    for oo=1:order
        
        if ieqtn==1
            
            perms{oo}=cell2mat(utils.gridfuncs.mypermutation(1:oo));
            
            initialize_output_order();
            
        end
        
        if verbose
            
            tic
            
        end
        
        ncols1=numel(myeqtn);
        
        myeqtn1=splanar.empty(1,0);
        
        ndvec=0;
        
        for icol1=1:ncols1
            
            obj=myeqtn(icol1);
            
            oldwrt=obj.lineage;
            
            ncols2=max(1,numel(oldwrt));
            
            pointer=0;
            
            for icol2=1:ncols2
                
                pointer=pointer+1;
                
                if ncols2>1
                    
                    obj2=intercept_column(obj,pointer);
                    
                else
                    
                    obj2=obj;
                    
                end
                
                candidates=find(obj2.incidence);
                
                if isempty(candidates)
                    
                    continue
                    
                end
                
                if isempty(oldwrt)
                    
                    b=[];
                    
                else
                    
                    b=oldwrt{pointer}; % lineage of derivatives...
                    
                end
                
                % Only consider variables that are greater than or equal to
                % the current pointing variable...
                wrt=candidates;
                
                if ~isempty(b)
                    
                    wrt(wrt<b(end))=[];
                    
                end
                
                if isempty(wrt)
                    
                    continue
                    
                end
                
                d=diff(obj2,wrt,alien_list);
                
                store_derivative()
                
            end
            
        end
        
        myeqtn=myeqtn1;
        
    end
    
end

    function store_derivative()
        
        if verbose
            
            derivs(oo).timer(ieqtn)=derivs(oo).timer(ieqtn)+toc;
            
        end
        
        % do not store zero derivatives
        %------------------------------
        dfunc=d.func;
        
        if isnumeric(dfunc) && isscalar(dfunc) && dfunc==0
            
            return
            
        end
        
        nw=numel(wrt);
        
        lin=cell(1,nw);
        
        locs=cell(1,nw+1);
        
        locs{1}=ieqtn;
        
        for ii=1:nw
            
            batch=[b,wrt(ii)];
            
            lin{ii}=batch;
            
            % do pairing
            %-----------            
            batch=batch(perms{oo});
            
            locs_i=utils.kronecker.remove_duplicated_tensors(batch);
            
            derivs(oo).nnz_derivs=derivs(oo).nnz_derivs+size(locs_i,1);
            
            locs{ii+1}=utils.gridfuncs.locate_permutation(locs_i,nwrt);
            
        end
        
        d.lineage=lin;
        
        d.location=locs;
        
        ndvec=ndvec+1;
        
        myeqtn1(ndvec)=d;
        
        derivs(oo).derivatives(end+1)=d;
        
    end

    function initialize_output_order()
                           
            derivs(oo).size=[neqtns,nwrt^oo];
            
            derivs(oo).nwrt=nwrt;
            
            derivs(oo).order=oo;
            
            derivs(oo).nnz_derivs=0;
            
            derivs(oo).derivatives=splanar.empty(0,1);
            
            if verbose
                
                derivs(oo).timer=zeros(neqtns,1);
                
            else
                
                derivs(oo).timer=nan(neqtns,1);
                
            end
                    
    end

end
