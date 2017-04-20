function out=decomp_format_output(out,y0,doLog)

fields=fieldnames(out);

nvar=get(y0,'NumberOfVariables');

is_multiple=nvar>1 && ~any(cellfun(@isempty,y0.varnames));

if is_multiple
    
    series_names=get(y0,'varnames');
    
    series_names=regexprep(series_names,'\s+','_');
    
    bad=~strcmp(y0.varnames,series_names);
    
    if any(bad)
        
        disp(y0.varnames(bad))
        
        warning('The variable names above were updated')
        
        y0.varnames(bad)=series_names(bad);
        
    end
    
    tmp=struct();
    
    y0p=pages2struct(y0);
    
    for sn=1:numel(series_names)
        
        cname=series_names{sn};
        
        tmp.(cname)=one_series(sn);
        
        tmp.(cname).y=y0p.(cname);
        
    end
    
    out=tmp;
    
else
    
    for kk=1:numel(fields)
        
        ff=fields{kk};
        
        dd=out.(ff);
        
        if doLog
            
            dd=exp(dd);
            
        end
        
        out.(ff)=reset_data(y0,dd);
        
    end
    
    out.y=y0;
    
end


    function thisSeries=one_series(position)
        
        thisSeries=struct();
        
        for jj=1:numel(fields)
            
            ff=fields{jj};
            
            dd_=out.(ff)(:,position);
            
            if doLog
                
                dd_=exp(dd_);
                
            end
            
            thisSeries.(ff)=reset_data(y0,dd_);
            
        end
        
    end

end