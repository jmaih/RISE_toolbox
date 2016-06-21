function [d,drop,start,summary]=process_draws(draws,drop)

if nargin<2
    
    drop=[];
    
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
    
    W=W.mat;
    
    if ~isempty(W)
        
        W=strrep(W,'.mat','');
        
        N=numel(W);
        
        re_order_names()
        
        d=cell(1,N);
        
        discard=false(1,N);
        
        for imat=1:N
            
            this_matrix=W{imat};
            
            tmp=load([draws,filesep,this_matrix]);
            
            if ~(isfield(tmp,'pop') && isfield(tmp.pop,'x'))
                
                discard(imat)=true;
                
                continue
                
            end
            
            d{imat}=tmp.pop;
            
            if isfield(tmp,'SIG')
                
                SIG=tmp.SIG;
                
            end
            
            if isfield(tmp,'c')
                
                c=tmp.c;
                
            end
            
        end
        
        d=d(~discard);
        
        d=[d{:}];
        
    end
    
else
    
    d=draws;
    
end

[nchains,npop]=size(d);

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
    'last_cScale',c);

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

    function re_order_names()
        
        Wbar=regexprep(W,'\w+_(\d+)','$1');
        
        for ii=1:N
            
            if isempty(Wbar{ii})
                
                Wbar{ii}=inf;
                
            else
                
                Wbar{ii}=str2double(Wbar{ii});
                
            end
            
        end
        
        Wbar=cell2mat(Wbar);
        
        [Wbar_ordered,tags]=sort(Wbar);
        
        last_saved_index=Wbar_ordered(end);
        
        W=W(tags);
        
    end

end