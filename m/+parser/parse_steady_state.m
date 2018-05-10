function [static,SteadyStateModel_block,auxiliary_steady_state_equations,...
    dictionary,blocks]=parse_steady_state(dictionary,blocks)

static=struct();

current_block_id=find(strcmp('steady_state_model',{blocks.name}));

if dictionary.parse_debug

    profile off
    
    profile on

else
    
    tic

end

[SteadyStateModel_block,dictionary,static]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'steady_state_model',static);

if dictionary.parse_debug

    profile off
    
    profile viewer
    
    keyboard

else
    
    disp([mfilename,':: Steady State Model parsing . ',sprintf('%0.4f',toc),' seconds'])

end

% get list of endogenous defined in the steady state: do this
% everytime there is an auxiliary variable
%----------------------------------------------------
nsstate=size(SteadyStateModel_block,1);

sstate_model_aux_vars=cell(1,nsstate);

iter=0;

for irow_=1:nsstate

    vname=SteadyStateModel_block{irow_,1}{1,1};
    
    if any(strcmp(vname,{dictionary.endogenous.name}))
    
        iter=iter+1;
        
        sstate_model_aux_vars{iter}=vname;
    
    end
    
end

dictionary.auxiliary_variables.ssmodel_solved=sstate_model_aux_vars(1:iter);

% remove item from block
blocks(current_block_id)=[];

auxiliary_steady_state_equations=dictionary.auxiliary_equations;

for irow_=1:size(auxiliary_steady_state_equations,1)

    tmp_=auxiliary_steady_state_equations{irow_,1}(1,:);
    
    tmp_{2}=strrep(tmp_{2},'-','=');
    
    auxiliary_steady_state_equations{irow_,1}(1,:)=tmp_;

end

end
