function H=set(H,varargin)

n=length(varargin);

if n==0
    
    return
    
end

if n==1
    
    arg=varargin{1};
    % struct style or query
    if isstruct(arg)
        % redo it
        do_struct2cell()
        
        return
    
    elseif ischar(arg)
        % query
        error('queries not implemented yet')
    
    else
        
        error('wrong input argument')
    
    end
    
elseif iscell(varargin{1})
    
    if ~iscell(varargin{2})|| n~=2
        
        error('wrong input argument')
    
    end
    
    arg1=varargin{1};  arg2=varargin{2};
    
    do_cell2cell()
    
    return
    
end

if rem(n,2)
    
    error('arguments must come in pairs')
    
end

H=setOptions(H,varargin{:});

    function do_cell2cell()
        
        [marg1,narg1]=size(arg1);
        
        if marg1~=1
            
            error('wrong input argument')
            
        end
        
        [marg2,narg2]=size(arg2);
        
        if narg1~=narg2
            
            error('wrong input argument')
            
        end
        
        mh=length(H);
        
        if marg2==1 && mh>1
            
            arg2=arg2(ones(1,mh),:);
            
            marg2=mh;
            
        end
        
        if marg2~=mh
            
            error('wrong input argument')
            
        end
        
        ff=[arg1;arg1];
        
        for ih=1:mh
            
            fi=ff;
            
            for icol=1:narg1
                
                fi{2,ih}=arg2{ih,icol};
                
            end
            
            fi=fi(:).';
            
            H(ih)=set(H(ih),fi{:});
            
        end
        
    end

    function do_struct2cell()
        
        ff=fieldnames(arg);
        
        ff=ff(:).';
        
        ff=[ff;ff];
        
        for ii=1:size(ff,2)
            
            ff{2,ii}=arg.(ff{1,ii});
            
        end
        
        ff=ff(:).';
        
        H=set(H,ff{:});
        
    end

end