function varargout=filter(self,param)

if nargin<2,param=[]; end

if isempty(param),param=self.estim_.estim_param; end

nout=nargout;

outCell=cell(1,nout);

%------------------------------------
XX=self.estim_.X;
                
YY=self.estim_.Y;

M=vartools.estim2states(param,...
    self.estim_.links.theMap,...
    self.mapping.nparams,...
    self.mapping.nregimes);

[outCell{1:nout}]=vartools.likelihood(M,self.mapping,...
    YY,XX,self.is_time_varying_trans_prob,...
    self.markov_chains);
%------------------------------------

% We need to shift the date to account for the fact that the VAR starts
% after its lags...
start_date=self.estim_.date_range(1)+self.nlags;

for ix=1:nout
    
    if isstruct(outCell{ix})
        
        outCell{ix}=reprocess_structure(outCell{ix});
        
    end
    
end

varargout=outCell;

    function sout=reprocess_structure(s)
        
        isshock=false;
        
        ff=fieldnames(s);
        
        sout=struct();
        
        for ii=1:numel(ff)
            
            sfi=s.(ff{ii});
            
            if self.is_panel
                
                ng=self.ng;
                
                isshock=isstruct(sfi) && ...
                    all(strncmp(fieldnames(sfi),'shock_',5));
                
                if isshock
                    
                    sfi=squash_shocks(sfi,ng);
                    
                end
                
                for g=1:ng
                    
                    sout.(self.members{g}).(ff{ii})=timeserize(sfi,g);
                    
                end
                
            else
                
                sout.(ff{ii})=timeserize(sfi);
                
            end
            
        end
        
        function sout=timeserize(s,g)
            
            if nargin<2
                
                g=[];
                
            end
            
            if isstruct(s)
                
                sout=struct();
                
                fff=fieldnames(s);
                
                for jj=1:numel(fff)
                    
                    sout.(fff{jj})=timeserize(s.(fff{jj}),g);
                    
                end
                
            else
                
                if ~isempty(g) && isshock
                    % chop
                    len=size(s,2);
                    
                    h=floor(len/ng);
                    
                    r=len-h*ng;
                    
                    s0=[];
                    
                    if r
                        
                        s0=s(:,1,:);
                        
                        s=s(:,2:end,:);
                        
                    end
                    
                    s=cat(2,s0,s(:,(g-1)*h+1:g*h,:));
                    
                end
                
                % would be nice to have one prototype with just
                % a start date...
                sout=ts(start_date,permute(s,[2,3,1]));
                
            end
            
        end
        
    end

end

function sfi=squash_shocks(sfi,ng)

if ng==1
    
    return
    
end

ffi=fieldnames(sfi);

if isstruct(sfi.(ffi{1}))
    
    return
    
end

nvars=numel(ffi);

nshocks=nvars/ng;

tmp=sfi.(ffi{1})(ones(nvars,1),:);

for ii=2:nvars
    
    tmp(ii,:)=sfi.(ffi{ii});
    
end

sfi=struct();

for ishock=1:nshocks
    
    v=tmp((ishock-1)*ng+1:ishock*ng,:,:);
    
   sfi.(ffi{ishock})=elongate(v);
    
end

    function e=elongate(v)
        
        e=cell(1,ng);
        
        for g=1:ng
            
            e{g}=v(g,:,:);
            
        end
        
        e=cat(2,e{:});
        
    end

end

%{
function sfi=squash_shocks(sfi,ng)

if ng==1
    
    return
    
end

ffi=fieldnames(sfi);

if isstruct(sfi.(ffi{1}))
    
    return
    
end

n=numel(ffi);

tmp=sfi.(ffi{1})(ones(n,1),:);

for ii=2:n
    
    tmp(ii,:)=sfi.(ffi{ii});
    
end

T=size(tmp,2);

nshocks=n/ng;

testla=reshape(tmp.',[T,ng,nshocks]);

ffi=ffi(1:nshocks);

sfi=struct();

for ii=1:nshocks
    
    sfi.(ffi{ii})=testla(:,:,ii).';
end

end

%}