function act=check_property(obj,propname,propval)

action=obj.spec_checker.(propname).action;

if nargout
    
    act=action;
    
end

if ~obj.spec_checker.(propname).check(propval)
    
    error(obj.spec_checker.(propname).errmsg)
    
end

end
