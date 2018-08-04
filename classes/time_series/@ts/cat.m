function db=cat(dim,varargin)
% cat concatenates time series along the specified dimension
%
% ::
%
%
%   db=cat(1,db1,db2,...,dbn)
%   db=cat(2,db1,db2,...,dbn)
%   db=cat(3,db1,db2,...,dbn)
%
% Args:
%
%    - **dim** [1|2|3] : dimension along which concatenation is done
%
%    - **db1, db2,...,dbn** [ts] : time series
%
% Returns:
%    :
%
%    - **db** [ts] : time series with concatenated series
%
% Note:
%
%    - all times series must be of the same frequency
%
%    - Concatenation along the second dimension requires that variables have
%      the same number of columns if no names are specified
%
%    - if names are specified in the first time series, then names should be
%      specified in all of the others as well.
%
%    - empty time series are discarded but there should be at least one
%      non-empty time series
%
% Example:
%
%    See also:


n=length(varargin);
dn1=cell(1,n);
nvar1=nan(1,n);
npages1=nan(1,n);
freq1=nan(1,n);
dn_max=-inf;
dn_min=inf;
discard=false(1,n);
for idb=1:n
    db_i=varargin{idb};
    if isempty(db_i)
        discard(idb)=true;
        continue
    end
    [dn1{idb},dn1_max,dn1_min,nvar1(idb),npages1(idb),freq1(idb)]=decompose_series(varargin{idb});
    dn_max=max(dn_max,dn1_max);
    dn_min=min(dn_min,dn1_min);
end
dn1=dn1(~discard);
nvar1=nvar1(~discard);
npages1=npages1(~discard);
freq1=freq1(~discard);
varargin=varargin(~discard);
n=sum(~discard);
if n==1
    db=varargin{1};
    return
end

if isempty(dn1)
    error('no valid time series to concatenate')
end
if ~all(freq1==freq1(1))
    error('data should have same frequency')
end
dn=dn_min:dn_max;
nobs=numel(dn);

no_names=isempty(varargin{1}.varnames{1});
switch dim
    case {1,3}
        if dim==1
            npages=max(npages1);
        else
            npages=sum(npages1);
        end
        [varnames,nvar]=set_variables_names();
        datta=nan(nobs,nvar,npages);
        offset_pages=0;
        for idb=1:n
            rows1=ts.set_locations(dn1{idb},dn);
            pp=offset_pages+(1:npages1(idb));
            if no_names
                datta(rows1,:,pp)=varargin{idb}.data;
            else
                varloc=locate_variables(varargin{idb}.varnames,varnames);
                datta(rows1,varloc,pp)=varargin{idb}.data;
            end
            if dim==3
                offset_pages=offset_pages+npages1(idb);
            end
        end
    case 2
        npages=max(npages1);
        % we do not insist on the names of the variables...
        varnames={};
        nvar=sum(nvar1);
        datta=nan(nobs,nvar,npages);
        offset_cols=0;
        for idb=1:n
            rows1=ts.set_locations(dn1{idb},dn);
            varnames=[varnames(:);varargin{idb}.varnames(:)];
            datta(rows1,offset_cols+(1:nvar1(idb)),1:npages1(idb))=varargin{idb}.data;
            offset_cols=offset_cols+nvar1(idb);
        end
    otherwise
        error('concatenation of time series allowed only for dimensions 1, 2 or 3')
end
db=ts(dn,datta,varnames);

    function [varnames,nvar]=set_variables_names()
        for idb_=1:n
            flag=isempty(varargin{idb_}.varnames{1});
            if no_names
                % no names are given. Then all guys should have the same number
                % of columns
                if varargin{idb_}.NumberOfVariables~=nvar1(1)
                    error('databases should have the same number of columns')
                end
                if idb_==1
                    varnames=varargin{1}.varnames;
                end
            else
                % all guys should have the names if the names are given. If they
                % don't all have the same names, the columns will be expanded
                if flag
                    error('at least one database has names and so all should have names')
                end
                if idb_==1
                    varnames=varargin{idb_}.varnames;
                else
                    varnames=union(varnames,varargin{idb_}.varnames);
                end
            end
        end
        nvar=numel(varnames);
    end
end