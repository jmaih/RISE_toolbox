function m=riff_erize(m)

m.routines.probs_times_dynamic=riff_erize_engine(m.routines.probs_times_dynamic);

for ii=1:numel(m.routines.probs_times_dynamic_derivatives)
    
m.routines.probs_times_dynamic_derivatives(ii).functions=...
    riff_erize_engine(m.routines.probs_times_dynamic_derivatives(ii).functions);

end

end

function out=riff_erize_engine(f)

mycells=cell(0,1); cellCount=0;

s0s1={}; ns0s1=[]; cc=[]; at_argins=[];

s0s1patt='\(?\(?\<s0\)?\s*==\s*\(?\d+\)?\s*&\s*\(?s1\)?\s*==\s*\(?\d+\)?\)?,?';

ifelseif='if_elseif';

cellstyle=iscell(f);


if cellstyle
    
    out=f;
    
    for ii=1:numel(f)
        
        out{ii}=main_engine(f{ii});
        
    end
    
else
    
    out={main_engine(f)};
    
end

memo=memo_(at_argins,s0s1,mycells);

out=functionalize(memo,out);

if ~cellstyle
    
    out=out{1};
    
end

    function out=main_engine(f)
        
        tmp=func2str(f);
        
        bingo=strfind(tmp,ifelseif);
        
        if isempty(bingo)
            
            out=f;
            
            return
            
        end
        
        first_closing=find(tmp==')',1,'first');
        
        at_argins=tmp(1:first_closing);
        
        parenth=sort([find(tmp=='('),find(tmp==')')]);
        
        for jj=numel(bingo):-1:1
            
            start=bingo(jj);
            
            pstart=find(parenth>start,1,'first');
            
            % plp=pstart;
            nopen=1;
            
            while nopen
                
                pstart=pstart+1;
                
                pos=parenth(pstart);
                
                if strcmp(tmp(pos),'(')
                    
                    nopen=nopen+1;
                    
                else
                    
                    nopen=nopen-1;
                    
                end
                
            end
            
            prp=pos;
            
            left=tmp(1:start-1);
            
            middle=tmp(start:prp);
            
            right=tmp(prp+1:end);
            
            tmp=[left,reprocess(middle),right];
            
            parenth=parenth(1:pstart-1);
            
        end
        
        out=tmp;
        
    end

    function out=reprocess(midl)
        
        midl=midl(length(ifelseif)+2:end-1);
        
        if isempty(s0s1)
            
            s0s1=regexp(midl,s0s1patt,'match');
            
            ns0s1=numel(s0s1);
            
            cc=cell(1,ns0s1);
                        
        end
        
        a=strfind(midl,s0s1{1});
        
        for iii=1:ns0s1
            
            al=length(s0s1{iii});
            
            if iii==ns0s1
                
                fn=midl(a+al:end);
                
            else
                
                b=strfind(midl,s0s1{iii+1});
                
                fn=midl(a+al:b-1);
                
                a=b;
                
            end
            
            cc{iii}=str2func([at_argins,fn]);
            
        end
        
        if isempty(mycells)
            
            mycells=cc;
            
        else
            
            mycells=[mycells;cc];
            
        end
        
        cellCount=cellCount+1;
        
        out=['memo(',int2str(cellCount),',',at_argins(3:end)];
        
    end

end

function out=functionalize(memo,out) %#ok<INUSL>

for ii=1:numel(out)
    
    if ischar(out{ii})
        
        out{ii}=eval(out{ii});% out{ii}=str2func(out{ii});
        
    end
    
end
        
end

function out=memo_(at_argins,s0s1,mycells)

mycond=str2func([at_argins,'[',cell2mat(s0s1),']']);

out=@thememoire;

    function out=thememoire(ping,varargin)
        
        out=riff(mycond,mycells(ping,:),varargin{:});
        
    end

end