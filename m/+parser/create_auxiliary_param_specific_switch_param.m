function [ssname,partname]=create_auxiliary_param_specific_switch_param(pname)
% create auxiliary parameters to be used in the perturbation approach by
% Foerster, Rubio-Ramirez, Waggoner and Zha(2016)

if nargin==0
    
    pname='';
    
end

if iscellstr(pname)
    
    n=numel(pname);
    
    ssname=cell(1,n);
    
    partname=ssname;
    
    for ii=1:n
        
        [ssname{ii},partname{ii}]=parser.create_auxiliary_param_specific_switch_param(pname{ii});
        
    end
    
else
    
    ssname=['ss_',pname];
    
    partname=['part_',pname];
    
end

end