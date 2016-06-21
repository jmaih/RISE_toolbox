function hdl=gelman_plot(obj,pname)

hdl0=plot(obj.psrf.(pname).time,obj.psrf.(pname).psrf);

title(sprintf('%s(%0.4f)',pname,obj.psrf.(pname).psrf(end)))

axis tight

if nargout
    
    hdl=hdl0;
    
end

end
