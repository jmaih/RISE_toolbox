function [xpand,keep,iter]=tensor_locator(nv,o,direction)

if nargin<3
    
    direction='down';
    
end

if ~ismember(direction,{'up','down'})
    
    error('third argument should be "up" or "down"')
    
end

if strcmp(direction,'up')
    
    direction_up=true;
    
else
    
    direction_up=false;
    
end


res=cell(1,o);

xpand=cell(1,o);

keep=cell(1,o);

npv=nv*ones(1,o);

for io=1:o
    
    res{io}=cell2mat(utils.gridfuncs.mypermutation(1:io));
    
    keep{io}=reshape(false(1,nv^io),[1,npv(1:io)]);
    
    xpand{io}=reshape(zeros(1,nv^io),[1,npv(1:io)]);
    
end

iter=zeros(1,o);

recurse()

    function recurse(bprev)
        
        if nargin==0
            
            bprev=[];
            
            if direction_up
                
                iprev=1;
                
            else
                
                iprev=nv;
                
            end
            
        else
            
            iprev=bprev(end);
            
        end
        
        curr_o=numel(bprev)+1;
        
        if direction_up
            
            range=iprev:nv;
            
        else
            
            range=1:iprev;
            
        end
        
        for icurr=range
            
            bcurr=[bprev,icurr];
            
            iter(curr_o)=iter(curr_o)+1;
            
            apply_position(bcurr,curr_o)
            
            if curr_o<o
                
                recurse(bcurr)
                
            end
            
        end
        
    end

    function apply_position(b,o)
        
        bx=b(res{o});
        
        pos=utils.gridfuncs.locate_permutation(bx,nv);
        
        keep{o}(pos(1))=true;
        
        %         pos=unique(pos);
        %         overwriting the same entry several times, surprisingly,
        %         turns out to be faster than finding the unique locations
        %         and writing to these.
        
        xpand{o}(pos)=iter(o);
        
    end

end
