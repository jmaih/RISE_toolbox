function plot_sims(m,vnames,range,varargin)

n=numel(vnames);

nr=floor(sqrt(n));

nc=ceil(sqrt(n));

odd=true;

while nr*nc<n
    if odd
        nr=nr+1;
    else
        nc=nc+1;
    end
    odd=~odd;
end

tex=get(m,'tex(long)');

figure;

islegend=false;

ndb=length(varargin);

leg=cell(1,ndb);

for ii=1:n
    
    subplot(nr,nc,ii)
    
    for jj=1:ndb
        
        tmp=varargin{jj}.(vnames{ii});
        
        if jj==1
            
            d=tmp;
            
        else
            
            d=[d,tmp]; %#ok<AGROW>
            
        end
        
        if ii==1
            
            if isfield(varargin{jj},'legend')
                
                leg{jj}=varargin{jj}.legend;
                
            end
            
            if jj==ndb
                
                islegend=ndb>1 && ~any(cellfun(@isempty,leg));
                
            end
            
        end
        
    end
    
    plot(range,d,'linewidth',2)
    
    title(tex.(vnames{ii}))
    
    if islegend && ii==1
        
        legend(leg,'Location','southeast')
        
        legend('boxoff')
        
    end
    
end

