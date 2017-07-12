function [Model_block,dictionary,blocks]=parse_model(dictionary,blocks)

current_block_id=find(strcmp('model',{blocks.name}));

more_string='';

if dictionary.definitions_inserted

    more_string='(& possibly definitions insertions)';

end

if dictionary.parse_debug

    profile off

    profile on

else
    
    tic

end

[Model_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'model');

if dictionary.parse_debug

    profile off
    profile viewer
    keyboard

else
    
    disp([mfilename,':: Model parsing ',more_string,'. ',sprintf('%0.4f',toc),' seconds'])

end

if isempty(Model_block)

    error([mfilename,':: no model declared'])

end
% remove item from block
blocks(current_block_id)=[];

% make auxiliary equations out of the list of endogenous variables
[dictionary,Model_block]=...
    parser.create_auxiliary_equations(dictionary,Model_block);

% Then modify the system for hybrid expectations
%------------------------------------------------
[Model_block,dictionary]=parser.hybrid_expectator(Model_block,dictionary);

neqtns=sum(strcmp(Model_block(:,end-1),'normal'));

nendo=numel(dictionary.endogenous);

if nendo<neqtns

    error(['More equations (',int2str(neqtns),') than the number of ',...
        'endogenous variables (',int2str(nendo),') '])

end

dictionary.is_deficient_eqtns=neqtns<nendo;

end
