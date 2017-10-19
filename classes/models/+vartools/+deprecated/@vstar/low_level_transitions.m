function [G,retcode,obj]=get_transitions(obj)

if isempty(obj)
    
    G=struct();
    
    return
    
end

obj=solve(obj);

thresholds_data=obj.data(obj.variables_locations_in_data.thresh_id,:);

retcode=0;

nthresh_=size(thresholds_data,1);

for icol=1:nthresh_
    
    sij=thresholds_data(icol,:);
    
    g=obj.solution.thresholds{icol}.g;
    
    c=num2cell(obj.solution.thresholds{icol}.c);
    
    if obj.options.estim_scale_free
        
        sd=std(sij);
        
        howmany=numel(c);
        
        if howmany==1 && ...
                strcmp(func2str(obj.solution.thresholds{icol}.func),...
                'exponential')
            
            howmany=2;
            
        end
        
        g=g*sd^howmany;
        
    end
    
    [G_icol,retcode]=obj.thresholds(icol).func(sij,g,c{:});
    
    if icol==1
        
        G=G_icol(ones(nthresh_,1),:);
        
    end
    
    G(icol,:) = G_icol;
    
    if retcode
        
        break
        
    end
    
end

G=G(:,obj.nlags+1:end);

end
