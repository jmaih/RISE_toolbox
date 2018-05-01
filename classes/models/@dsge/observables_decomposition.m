function varargout=observables_decomposition(obj,select,xrange,varargin)
% OBSERVABLES_DECOMPOSITION - decomposes all variables of a DSGE model in
% terms of observables.
%
% ::
%
%
%   weights=observables_decomposition(obj,select,xrange,db)
%
%   [weights,dec1,...,decn]=observables_decomposition(obj,select,xrange,db1,...,dbn)
%
% Args:
%
%    - **obj** [rise|dsge]: scalar model object
%
%    - **select** [{[]}|'a'|'att'|'alpha'|'v'|'r'|'epsilon'|'eta']: type of
%    decomposition to perform
%      - **a** : filter
%      - **att** : update
%      - **alpha** : smooth
%      - **v** : forecast errors
%      - **r** : variables calculated during the smoothing process
%      - **epsilon** : measurement errors
%      - **eta** : structural shocks
%      - **[]** : all the above
%
%    - **xrange** [{[]}|vector|cell array]: start date and end date for the
%    decomposition
%
%    - **db** [ts]: database with the data to be used in the decomposition
%
% Returns:
%    :
%
%    - **weights** [struct]: weights for the different elements requested from
%    the variable **select**
%
%    - **dec1** [struct]: hyper structure containing the time series for the
%    various decomposition types and variables
%
% Note:
%
%    - if **select** is empty, all the decompositions are performed
%
%    - After doing the decomposition of different variables, one can take the
%    differences in the decomposition to perform, e.g., the decomposition of
%    forecast errors in terms of observables.
%
% Example:
%
%    See also: HISTORICAL_DECOMPOSITION

if isempty(obj)
        
    if nargout
        
        varargout={struct.empty};
                
    end
    
    return
    
end

do_demean=true;

nobj=numel(obj);

if nobj>1
    
    error('this function cannot be vectorized')
    
end

ndatasets=numel(varargin);

if ndatasets==0
    
    error('at least one dataset should be provided')
    
end

dataOut=cell(1,ndatasets);

if ~isempty(xrange)
    
    if isa(xrange,'double')
        
        xrange=num2cell(xrange);
    
    end
    
    if numel(xrange)<2
        
        error('third input for the range should have at least two elements')
    
    end
    
    obj=set(obj,'estim_start_date',xrange{1},'estim_end_date',xrange{end});
    
end

for ii=1:ndatasets
    
    [dataOut{ii},~,start_date,~,retcode]=data_prerequest(obj,varargin{ii});
    
    if retcode
        
        error(decipher(retcode))
        
    end
    
    dataOut{ii}=dataOut{ii}.y;
    
end

[obj,retcode]=solve(obj);

if retcode
    
    if obj(1).options.debug
        
        error(decipher(retcode))
        
    end
    
    return
    
end

[syst,retcode]=filter_initialization(obj);

if retcode
    
    error(decipher(retcode))
    
end

nout=nargout;

[varargout{1:nout}]=main_engine(syst,select,do_demean,obj.options.debug,dataOut{:});

obsnames=obj.observables.name;
vnames=struct();
vnames.a=obj.endogenous.name(obj.order_var);
vnames.att=vnames.a;
vnames.alpha=vnames.a;
vnames.r=vnames.a;
vnames.epsilon=obsnames;
vnames.eta=obj.exogenous.name;
vnames.v=obsnames;

proto=ts(start_date,zeros(1,obj.observables.number(1)));

for jj=2:nout
    
    varargout{jj}=time_serize(varargout{jj});
    
end

    function vout=time_serize(vin)
        
        fields=fieldnames(vin);
        
        vout=struct();
        
        for iii=1:numel(fields)
        
            ff=fields{iii};
            
            inames=vnames.(ff);
            
            for jjj=1:numel(inames)
                
                newdata=vin.(ff)(:,:,jjj);
            
                vout.(fields{iii}).(inames{jjj})=reset_data(proto,...
                    newdata,obsnames);
            
            end
            
        end
            
    end

end

function varargout=main_engine(syst,select,do_demean,debug,varargin)

ndatasets=length(varargin);

if nargout>ndatasets+1
    
    error('too many output arguments')
    
end

[T,R,Z,H,Q,sstate,init,growth]=dsge.state_space_wrapper(syst);

first_dataset=varargin{1};

[nvobs,nobs]=size(first_dataset);

for id=2:ndatasets
    
    if ~isequal(size(varargin{id}),[nvobs,nobs])
        
        error('all datasets should have the same size')
        
    end
    
end

ybar=sstate(Z);

dy=ybar-sstate(Z);

m=size(T,1);

ca=(eye(m)-T)*sstate;

if any(growth~=0)
    
    ca = ca+growth;
    
end

if do_demean
    
    first_dataset=bsxfun(@minus,first_dataset,sstate(Z));
    
    ca=zeros(size(ca));
    
    dy=zeros(size(dy));
    
    init.a=init.a-sstate;
    
end

ff=obsw.missing_observations_kalman_filter(first_dataset,T,R,Z,H,Q,init,dy,ca);

if ~debug
    
    first_dataset=[];
    
end

w=obsw.observation_weights(Z,T,H,Q,R,ff,select,first_dataset);

% update selection in case it was empty above
%--------------------------------------------
select=fieldnames(w);

varargout=cell(1,ndatasets+1);

varargout{1}=w;

m=size(T,1);

nv=struct('a',m,'att',m,'alpha',m,'v',nvobs,'r',m,'epsilon',nvobs,'eta',size(R,2));

for id=1:ndatasets
    
    datai=varargin{id};
    
    missing=isnan(datai);
    
    if id==1
        
        missing1=missing;
        
    else
        
        if any(vec(missing-missing1)~=0)
            
            error('datasets should have the same missing structure')
            
        end
        
    end
    
    if do_demean
        
        datai=bsxfun(@minus,datai,sstate(Z));
        
    end
    
    varargout{1+id}=do_one(datai);
    
end

    function out=do_one(y)
        
        % ensure that the multiplication works
        %--------------------------------------
        y(isnan(y))=0;
        
        out=struct();
        
        for ii=1:numel(select)
            
            field=select{ii};
            
            out.(field)=do_one_field(field);
            
        end
        
        function VV=do_one_field(field)
            
            way=bsxfun(@times,w.(field),y(:).');
            
            [r,c]=size(way);
            
            % aggregate contributions of each variable across time
            n=c/nvobs;
            
            V=0;
            
            for jj=1:n
                
                cols=(jj-1)*nvobs+1:jj*nvobs;
                
                V=V+way(:,cols);
                
            end
            
            nvf=nv.(field);
            
            n=r/nvf;
            
            VV=zeros(n,nvobs,nvf);
            
            % separate variables
            %-------------------
            for iv=1:nvf
                
                stretch=iv:nvf:r;
                
                VV(:,:,iv)=V(stretch,:);
                
            end
            
        end
        
    end

end