function [d,drop,start,summary]=process_draws(draws,drop,start_from,trimming)

if nargin<4
    
    trimming=[];
    
    if nargin<3
        
        start_from=[];
        
        if nargin<2
            
            drop=[];
            
        end
        
    end
    
end

if isempty(trimming)
    
    trimming = 1;
    
end

if isempty(start_from)
    
    start_from=1;
    
end

if isempty(drop)
    
    drop=50/100;
    
else
    
    if drop<0 || drop >=1
        
        error('drop should be in [0,1)')
        
    end
    
end

is_saved_to_disk=ischar(draws);

last_saved_index=0;

SIG=[];

d=[];

c=[];

if is_saved_to_disk
    
    W = what(draws);
    
    if ~isempty(W)
        
        W=W.mat;
        
        W=strrep(W,'.mat','');
        
        N=re_order_names();
        
        d=cell(1,N);
        
        iter=0;
        
        offset=0;
        
        for imat=1:N
            
            this_matrix=W{imat};
            
            tmp=load([draws,filesep,this_matrix]);
            
            if ~(isfield(tmp,'pop') && isfield(tmp.pop,'x'))
                
                continue
                
            end
            
            iter=iter+1;
            
            nc=size(tmp.pop,2);
            
            if iter==1
                
                
                d=tmp.pop(:,ones(1,1e+6));
                
            end
            
            if offset+nc>size(d,2)
                
                d=[d,tmp.pop(:,ones(1,1e+5))]; %#ok<AGROW>
                
            end
            
            d(:,offset+(1:nc))=tmp.pop;
            
            offset=offset+nc;
            
            if isfield(tmp,'SIG')
                
                SIG=tmp.SIG;
                
            end
            
            if isfield(tmp,'c')
                
                c=tmp.c;
                
            end
            
        end
        
        d=d(:,1:offset);
        
    end
    
else
    
    d=draws;
    
end

if isa(trimming,'function_handle')
    
    d=trimming(d);
    
else
    
    d=d(:,1:trimming:end);

end

[nchains,npop]=size(d);

if npop==0
    
    nchains=0;
    
end

discard=round(drop*npop);

start=discard+1;

% store the best before proceeding
%---------------------------------
best=load_best();

d=d(:,start:end);

dlast=[];

% estimated covariances
%-----------------------

for ichain=1:nchains
    
    theCov=cov([d(ichain,:).x].');
    
    if ichain==1
        
        eSIG=theCov;
        
        eSIG=eSIG(:,:,ones(nchains,1));
        
    else
        
        eSIG(:,:,ichain)=theCov;
        
    end
    
    if ichain==nchains
        
        dlast=d(:,end);
        
    end
    
end

if nchains==0
    
    eSIG=[];
    
end

summary=struct('nchains',nchains,'npop',npop,...
    'last_saved_index',last_saved_index,...
    'best_of_the_best',best,'last',dlast,...
    'last_cov',SIG,'estimated_cov',eSIG,...
    'last_cScale',c,...
    'trimming',trimming);

    function best=load_best()
        
        best=[];
        
        if isfield(d,'f')
            
            best=1;
            
            bestf=d(best).f;
            
            for iii=2:numel(d)
                
                if d(iii).f<bestf
                    
                    best=iii;
                    
                    bestf=d(best).f;
                    
                end
                
            end
            
            best=d(best);
            
        end
        
    end

    function N=re_order_names()
        
        if isempty(W)
            
            N = 0;
            
            return
            
        end
        
        Wbar=regexprep(W,'\w+_(\d+)','$1');
        
        for ii=1:numel(Wbar)
            
            if isempty(Wbar{ii})
                
                Wbar{ii}=inf;
                
            else
                
                Wbar{ii}=str2double(Wbar{ii});
                
            end
            
        end
        
        Wbar=cell2mat(Wbar);
        
        [Wbar_ordered,tags]=sort(Wbar);
        
        last_saved_index=Wbar_ordered(end);
        
        W=W(tags(start_from:end));
        
        N=numel(W);
        
    end

end