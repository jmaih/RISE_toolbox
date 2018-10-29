function [blocks,logvar_names]=file2blocks_log_variables(blocks)
% INTERNAL FUNCTION
%

endogBlock=strcmp('endogenous',{blocks.name});
endovar_names={blocks(endogBlock).listing.name};

levelVarBlock=strcmp('level_variables',{blocks.name});
levelvar_names={blocks(levelVarBlock).listing.name};

logvarBlock=strcmp('log_vars',{blocks.name});
logvar_names=load_empty_log_variables();

if ~isempty(levelvar_names)
    
    if ~isempty(logvar_names)
        
        error('log_variables and level_variables cannot be declared in the same file')
        
    end
    
    store_locs=locator(levelvar_names);
    
    blocks(levelVarBlock).listing=blocks(logvarBlock).listing;
    
    % now swap and destroy
    blocks(logvarBlock).listing=blocks(endogBlock).listing;
    
    if blocks(levelVarBlock).is_all_but
        
        blocks(logvarBlock).listing=blocks(logvarBlock).listing(store_locs);
        
        blocks(levelVarBlock).is_all_but=false;
        
    else
        
        blocks(logvarBlock).listing(store_locs)=[];
        
    end
    
    logvar_names={blocks(logvarBlock).listing.name};
    
end

% check that all log_vars are endogenous
%---------------------------------------
if ~isempty(logvar_names)
    
    store_locs=locator(logvar_names);
    
    if blocks(logvarBlock).is_all_but
        
        blocks(logvarBlock).listing=blocks(endogBlock).listing;
        
        blocks(logvarBlock).listing(store_locs)=[];
        
        logvar_names={blocks(logvarBlock).listing.name};
        
    end
    
end

    function logvar_names=load_empty_log_variables()
        
        logvar_names={blocks(logvarBlock).listing.name};
        
        if isempty(logvar_names) && blocks(logvarBlock).is_all_but
            
            blocks(logvarBlock).listing=blocks(endogBlock).listing;
            
            blocks(logvarBlock).is_all_but=false;
            
            logvar_names={blocks(logvarBlock).listing.name};
            
        end
        
    end

    function store_locs=locator(vList)
        
        store_locs=locate_variables(vList,endovar_names,true);
        
        locs=find(isnan(store_locs));
        
        if ~isempty(locs)
            
            bad_vars=vList(locs);
            
            disp(bad_vars(:)')
            
            error('The LEVEL variables above have not been found in the list of endogenous variables')
            
        end
        
    end

end
