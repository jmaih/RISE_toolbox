function [d,fd,c,SIG,accept_ratio,last_saved_index,dlast,best]=reload_draws(simfold,subset)
% subset is a 1 x 2 cell array in which
% - the first cell contains the columns to retain: can be empty, defaults
% to all
% - the second column contains the chains to retain: can be empty, defaults
% to all
if nargin<2
    
    subset=[];
    
    if nargin<1
        
        error('the function should have at least one input argument')
        
    end
    
end

last_saved_index=0;

is_saved_to_disk=ischar(simfold);

if is_saved_to_disk
    
    preselected_chains=[];
    
    if ~isempty(subset) && numel(subset)==2
        
        preselected_chains=subset{2};
        
    end
    
    [d,fd,c,SIG,accept_ratio,dlast,best]=load_from_disk();
    
else
    
    [d,fd,c,SIG,accept_ratio,dlast,best]=load_from_self();
    
end

if isempty(d)
    
    return
    
end

[~,ncols,nchains]=size(d);

set_subset(ncols,nchains)

d=d(:,subset{1},subset{2});

fd=fd(1,subset{1},subset{2});

c=c(subset{2});

SIG=SIG(subset{2});

accept_ratio=accept_ratio(subset{2});

if iscell(last_saved_index)
    % only valid for saved to disk
    last_saved_index=last_saved_index(subset{2});
    
end

dlast=dlast(subset{2});

best=best(subset{2});


    function [d,ff,c,SIG,accept_ratio,dlast,best]=load_from_self()
        
        if isstruct(simfold)
            
            simfold={simfold};
            
        end
        
        nchains=numel(simfold);
        
        c=cell(1,nchains); SIG=c; accept_ratio=c;  dlast=c;  best=c;
        
        for ic=1:nchains
            
            dii=simfold{ic};
            
            if ic==1
                
                d=[dii.pop.x];
                
                d=d(:,:,ones(1,nchains));
                
                ff=d(1,:,:);
                
            end
            
            d(:,:,ic)=[dii.pop.x];
            
            ff(1,:,ic)=[dii.pop.f];
            
            SIG{ic}=dii.SIG;
            
            c{ic}=dii.c;
            
            accept_ratio{ic}=dii.stats.accept_ratio;
            
            dlast{ic}=dii.pop(end).x;
            
            best{ic}=dii.best;
            
        end
        
    end


    function [d,ff,c,SIG,accept_ratio,dlast,best]=load_from_disk()
        
        a=what(simfold);
        
        a=strrep(a.mat,'.mat','');
        
        if isempty(a)
            
            d={}; c={}; SIG={}; accept_ratio={}; dlast={}; best={}; 
            
            return
            
        end
        
        proto=load([simfold,filesep,a{1}]);
        
        np=numel(proto.bestx);
        
        npop=numel(proto.pop);
        
        a=regexp(a,'(?<stud>mhDraws)_(?<chain>\d+)_(?<batch>\d+)','names');
        
        a=[a{:}];
        
        n=numel(a);
        
        for v0={'chain','batch'}
            
            v=v0{1};
            
            cc=cellfun(@(x)str2double(x),{a.(v)},'uniformOutput',false);
            
            [a(1:n).(v)]=deal(cc{:});
            
        end
        
        mhDraws=a(1).stud;
        
        chains=unique([a.chain]);
        
        nchains=numel(chains);
        
        nmat=max([a.batch]);
        
        d=zeros(np,npop*nmat,nchains);
        
        ff=zeros(1,npop*nmat,nchains);
        
        offset=zeros(1,nchains);
                
        c=cell(1,nchains); SIG=c; accept_ratio=c; 
        last_saved_index=c; dlast=c; best=c;
        
        if isempty(preselected_chains)
            
            preselected_chains=1:nchains;
            
        end
        
        for ic=1:nchains
            
            if ~any(preselected_chains-ic==0)
                % don't waste time loading chains that are not needed
                offset(ic)=inf;
                
                continue
                
            end
            
            athis=a([a.chain]==ic);
            
            nmat_ic=max([athis.batch]);
            
            for imat=1:nmat_ic
                
                file=sprintf('%s%s%s_%0.0d_%0.0d',simfold,filesep,mhDraws,ic,imat);
                
                dii=load(file);
                
                npi=numel(dii.pop);
                
                d(:,offset(ic)+(1:npi),ic)=[dii.pop.x];
                
                ff(:,offset(ic)+(1:npi),ic)=[dii.pop.f];
                
                offset(ic)=offset(ic)+npi;
                
                % recovering of simulations
                %---------------------------
                if imat==nmat_ic
                    
                    SIG{ic}=dii.SIG;
                    
                    c{ic}=dii.c;
                    
                    accept_ratio{ic}=dii.stats.accept_ratio;
                    
                    last_saved_index{ic}=imat;
                    
                    dlast{ic}=dii.pop(end).x;
                    
                    best{ic}=dii.best;
                    
                end
                
            end
            
        end
        
        d=d(:,1:min(offset),:);
        
        ff=ff(1,1:min(offset),:);
        
    end


    function set_subset(ncols,nchains)
        
        if isempty(subset)
            
            subset=cell(1,2);
            
        elseif iscell(subset)
            
            if numel(subset)==1
                
                subset=[subset,{[]}];
                
            end
            
        else
            
            error('subset must be a cell array')
            
        end
        
        binge=[ncols,nchains];
        
        for ii=1:2
            
            if isempty(subset{ii})
                
                subset{ii}=1:binge(ii);
                
            end
            
            subset{ii}=sort(subset{ii});
            
            if ii==1
                
                if numel(subset{ii})==2
                    
                    subset{ii}=subset{ii}(1):subset{ii}(2);
                    
                end
                
            end
            
            if subset{ii}(1)>binge(ii)||subset{ii}(end)>binge(ii)
                
                error('numbers in subsetting incompatible with available data')
                
            end
            
        end
        
    end

end
