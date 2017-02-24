function [checker,myoptions]=set_checker(checker,myoptions,new)

if isempty(checker),checker=struct(); end

if isempty(myoptions),myoptions=struct(); end

for ii=1:size(new,1)
    
    engine(new(ii,:))
    
end

    function engine(novel)
        
        prop=novel{1};
        
        action='';
        
        lp=find(prop=='(');
        
        if ~isempty(lp)
            
            action=prop(lp+1:end-1);
            
            prop=prop(1:lp-1);
            
        end
        
        checker.(prop).default=novel{2};
        
        myoptions.(prop)=checker.(prop).default;
        
        checker.(prop).check=novel{3};
        
        checker.(prop).errmsg=novel{4};
        
        checker.(prop).action=action;
        
    end


end