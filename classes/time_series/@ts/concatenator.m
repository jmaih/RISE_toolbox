function myirfs=concatenator(varargin)
% CONCATENATOR - concatenate outputs from multiple time series or
% structures of time series
%
% ::
%
%
% - myirfs=CONCATENATOR(db1,db2,...,dbn)
%
% Args:
%
%    - **dbi** [struct|ts] : time series or structure of time series
%
% Returns:
%    :
%
%    - **myirfs** [struct] : structure containing the concatenated time series
%
% Note:
%
% Example:
%
%    See also:

myirfs=struct();

n=length(varargin);

[fnames,is_ts]=gather_relevant_fields();

for ii=1:numel(fnames)
    
    pointer=fnames{ii};
    
    newbatch=cell(1,n);
    
    for jj=1:n
        
        newbatch{jj}=varargin{jj}.(pointer);
        
    end
    
    if is_ts
        
        myirfs.(pointer)=do_it(newbatch{:});
        
    else
        
        myirfs.(pointer)=ts.concatenator(newbatch{:});
        
    end
    
end

    function [fnames,is_ts]=gather_relevant_fields()
        
        for jjj=1:n
            
            if isa(varargin{jjj},'ts')
                
                varargin{jjj}=pages2struct(varargin{jjj});
                
            end
            
            if jjj==1
                
                fnames=fieldnames(varargin{1});
                
                is_ts=isa(varargin{1}.(fnames{1}),'ts');
                
                continue
                
            end
            
            fj=fieldnames(varargin{jjj});
            
            fnames_=intersect(fnames,fj);
            
            if isempty(fnames_)
                
                error('no common variables')
                
            end
            
            d=fnames-fnames_;
            
            if ~isempty(d)
                
                disp(d)
                
                warning(['the following fields do not belong to the intersection ',...
                    'and are therefore discarded'])
                
            end
            
            fnames=fnames_;
            
        end
        
    end

end

function myts=do_it(varargin)

nn=length(varargin);

vnames=varargin{1}.varnames;

nnames=numel(vnames);

if nnames==1
    
    [myts,dn,nobs,np]=dissect(varargin{1});
    
    myts=myts(:,ones(1,nn),:);
    
    for jjj=2:nn
        
        stretch_batch(jjj)
        
    end
    
    myts=reset_data(varargin{1},myts);
    
else
    
    batch=cell(1,nn);
    
    myts=struct();
    
    s=struct('type','()','subs','hallo');
    
    for jjj=1:nnames
        
        v=vnames{jjj};
        
        for iii=1:nn
            
            %batch{iii}=varargin{iii}(v);
            s.subs={v};
            batch{iii}=subsref(varargin{iii},s);
            
        end
        
        myts.(v)=do_it(batch{:});
        
    end
    
end

    function stretch_batch(id)
        
        [dataj,dnj,~,npj]=dissect(varargin{id});
        
        npdiff=npj-np;
        
        if npdiff>0
            
            myts=cat(3,myts,nan(nobs,nn,npdiff));
            
            np=npj;
            
        end
        
        dnstart=dn(1); while dnstart>dnj(1),dnstart=dnstart-1; end
        
        dnfinish=dn(end); while dnfinish<dnj(end),dnfinish=dnfinish+1; end
        
        start_diff=dn(1)-dnstart;
        
        if start_diff>0
            
            myts=cat(1,nan(start_diff,nn,np),myts);
            
        end
        
        finish_diff=dnfinish-dn(end);
        
        if finish_diff>0
            
            myts=cat(1,myts,nan(finish_diff,nn,np));
            
        end
        
        dn=dnstart:dnfinish;
        
        nobs=numel(dn);
        
        a=find(dn==dnj(1));
        
        b=find(dn==dnj(end));
                
        myts(a:b,id,1:npj)=dataj;
        
    end

    function [data,date_numbers,nobs,np]=dissect(tidsserie)
        
        data=tidsserie.data;
        
        date_numbers=tidsserie.date_numbers;
        
        nobs=tidsserie.NumberOfObservations;
        
        np=tidsserie.NumberOfPages;
        
    end

end
