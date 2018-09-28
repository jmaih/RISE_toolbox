function [shadow_model,sig]=add_auxiliary_switching_parameters(shadow_model,dicparams)
% replace all switching parameters with an auxiliary expression involving
% the perturbation parameter and the steady state for the switching
% parameter

sp=[dicparams.is_switching];

sig='';

if ~any(sp)
    
    return
    
end

pnames={dicparams.name};

sp=pnames(sp);

n=numel(sp);

[ssnames,partnames]=parser.create_auxiliary_param_specific_switch_param(sp);

others=parser.perturbation_control_param_names();

allswitch=[sp(:).',ssnames(:).',partnames(:).',others(:).'];

pos=locate_variables(allswitch,pnames);

panal=parser.parameters_analytical_form(pos);

sigpos=3*n+1;

iota1pos=sigpos+1;

sig=panal{sigpos};

iota1=panal{iota1pos};

switch_map=reshape(panal(1:n*3),[],3);

do_replace=@this_is_slick; %#ok<NASGU>

express=['\<',parser.cell2matize(panal(1:numel(sp))),'\>'];

express=regexprep(express,'\((\d+)\)','\\($1\\)');

shadow_model=regexprep(shadow_model,express,'${do_replace($1)}');

    function out=this_is_slick(v)
        
        loc=strcmp(switch_map(:,1),v);
        
        thetabar=switch_map{loc,2};
        
        pj=switch_map{loc,3};
        
        out=['frwzsp(',v,',',thetabar,',',pj,',',iota1,',',sig,')'];
        
    end


end