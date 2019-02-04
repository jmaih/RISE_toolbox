function [d,drop,start,summary]=process_draws(draws,whichChains,drop,...
    start,trimming)
% INTERNAL FUNCTION
%

if nargin<5
    
    trimming=[];
    
    if nargin<4
        
        start=[];
        
        if nargin<3
            
            drop=[];
            
            if nargin<2
                
                whichChains=[];
                
            end
            
        end
        
    end
    
end

if isempty(trimming),trimming = 1; end

if isempty(start)
    
    start=1;
    
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
        
        [N,chains,stud]=re_order_names();
        
        nchains=numel(N);
        
        if isempty(whichChains)||nchains==0
            
            whichChains=chains;
            
        elseif nchains && ~all(ismember(whichChains,chains))
            
            nchains=0;
            
            whichChains=1:nchains;
            
        end
        
        if nchains
            
            last_saved_index=last_saved_index(whichChains);
            
            nchains=numel(whichChains);
            
        end
        
        d=cell(1,nchains); c=cell(1,nchains); SIG=cell(1,nchains);
        accept_ratio=cell(1,nchains);
        
        iter=0;
        
        for ichain=whichChains
            
            iter=iter+1;
            
            [d{iter},c{iter},SIG{iter},accept_ratio{iter}]=load_chain(stud,...
                ichain,N(ichain));
            
        end
        
    end
    
else
    
    if isstruct(draws)
        
        draws={draws};
        
    end
    
    nchains=numel(draws);
    
    if isempty(whichChains)
        
        whichChains=1:nchains;
        
    end
    
    nchains=numel(whichChains);
    
    draws=draws(whichChains);
    
    d=cell(1,nchains); c=cell(1,nchains); SIG=cell(1,nchains);
    last_saved_index=cell(1,nchains);
    accept_ratio=cell(1,nchains);
    
    for ichain=1:nchains
        
        SIG{ichain}=draws{ichain}.SIG;
        
        c{ichain}=draws{ichain}.c;
        
        d{ichain}=draws{ichain}.pop;
        
        % index of mat files saved to disk
        %---------------------------------
        last_saved_index{ichain}=0;
        
        accept_ratio{ichain}=draws{ichain}.stats.accept_ratio;
        
    end
    
end

if isa(trimming,'function_handle')
    
    d=cellfun(trimming,d,'uniformOutput',false);
    
else
    
    d=cellfun(@(x)x(:,1:trimming:end),d,'uniformOutput',false);
    
end

% store the best before proceeding
%---------------------------------
best=load_best();

npop=cellfun(@(x)size(x,2),d,'uniformOutput',false);

dlast=cell(1,nchains);

eSIG=cell(1,nchains);

% estimated covariances
%-----------------------
for ichain=1:nchains
    
    discard=round(drop*npop{ichain});
    
    start=discard+1;
    
    if isempty(d{ichain})
        
        continue
        
    end
    
    d{ichain}=d{ichain}(:,start:end);
    
    %     theCov=cov([d{ichain}.x].');
    
    %     eSIG{ichain}=cov([d{ichain}.x].');
    
    dlast{ichain}=d{ichain}(:,end);
    
end

summary=struct('nchains',nchains,'npop',npop,...
    'last_saved_index',last_saved_index,...
    'best_of_the_best',best,'last',dlast,...
    'last_cov',SIG,'estimated_cov',eSIG,...
    'last_cScale',c,...
    'trimming',trimming,...
    'accept_ratio',accept_ratio);

    function best=load_best()
        
        best=cell(1,nchains);
        
        for jchain=1:nchains % if isfield(d,'f')
            
            if isempty(d{jchain})
                
                continue
                
            end
            
            [~,logic]=min([d{jchain}.f]);
            
            best{jchain}=d{jchain}(logic);
            
        end
        
    end

    function [N,chains,stud]=re_order_names()
        
        if isempty(W)
            
            N = [];
            
            stud='';
            
            chains=[];
            
            return
            
        end
        
        WS=regexp(W,'(?<stud>[^_]+)_(?<chain>\d+)_(?<matfile>\d+)','names');
        
        WS=[WS{:}];
        
        stud=WS(1).stud;
        
        specchains=cellfun(@(x)str2double(x),{WS.chain});
        
        matfiles=cellfun(@(x)str2double(x),{WS.matfile});
        
        chains=unique(specchains);
        
        nchain=numel(chains);
        
        N=zeros(1,nchain);
        
        for ichene=1:nchain
            
            current=specchains==ichene;
            
            if any(current)
                
                N(ichene)=max(matfiles(current));
            
            end
            
        end
        
        last_saved_index=num2cell(N);
        
    end

    function [d,c,SIG,accept_ratio]=load_chain(stud,chain,N)
        
        d=[]; c=[]; SIG=[];
        
        offset=0;
        
        accept_ratio=0;
        
        for jmat=1:N
            
            this_matrix=sprintf('%s_%d_%d',stud,chain,jmat);
            
            tmp=load([draws,filesep,this_matrix]);
            
            if ~(isfield(tmp,'pop') && isfield(tmp.pop,'x'))
                
                continue
                
            end
            
            nc=size(tmp.pop,2);
            
            if isempty(d)
                
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
            
            if jmat==N
                
                accept_ratio=tmp.stats.accept_ratio;
                
            end
            
        end
        
        d=d(:,1:offset);
        
    end

end