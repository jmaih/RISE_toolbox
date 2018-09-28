function out=if_elseif_to_if_elseiff(in)

if iscell(in)
    
    out=in;
    
    for ii=1:numel(in)
        
        out{ii}=parser.if_elseif_to_if_elseiff(out{ii});
        
    end
    
    return
    
end

fh=isa(in,'function_handle');

if fh
    
    fstr0=func2str(in);
    
else
    
    fstr0=in;
    
end

% remove calls to if_elseif
%--------------------------
[start,finish]=regexp(fstr0,'if_elseif','start','end');

if isempty(start)
    
    out=in;
    
    return
    
end

preamble=fstr0(1:start(1)-1);

first_close=find(preamble==')',1,'first');

args=preamble(3:first_close-1);

for ii=numel(start):-1:1
    
    process_one(start(ii),finish(ii));
    
end

out=fstr0;

if fh
    
    out=str2func(out);
    
end

    function process_one(start,last)
        
        op=1;
        
        iter=last+1;
        
        while op
            
            iter=iter+1;
            
            atom=fstr0(iter);
            
            if atom=='('
                
                op=op+1;
                
            elseif atom==')'
                
                op=op-1;
                
            end
            
        end
        
        ifelsiff_=fstr0(start:iter);
        
        rest=fstr0(iter+1:end);
        
        studs=regexp(ifelsiff_,'s0==\d+&s1==\d+,','split');
        % ignore if_elseif
        studs=studs(2:end);
        % remove the last element in all entries (commas and closing parenthesis)
        studs1=cellfun(@(x)x(1:end-1),studs,'uniformOutput',false);
        
        voila=regexp(ifelsiff_,'s0==\d+&s1==\d+','match');
        
        studs1=strcat(['@(',args,')'],studs1);
        
        key=[voila(:).';studs1(:).'];
        
        key=cell2mat(strcat(key(:).',','));
        
        newdeal=['if_elseiff({',args,'},',key(1:end-1),')'];
        
        fstr0=[fstr0(1:start-1),newdeal,rest];
        
    end


end