function [dictionary,blocks]=parse_exogenous_definitions(dictionary,blocks)
% this block needs a special treatment as the information provided here
% will only be used when building the dataset for estimation and/or during
% forecasting.
current_block_id=find(strcmp('exogenous_definition',{blocks.name}));
if dictionary.parse_debug
    profile off
    profile on
else
    tic
end
[ExogenousDefinition_block,dictionary]=parser.capture_equations(dictionary,blocks(current_block_id).listing,'exogenous_definition');
if dictionary.parse_debug
    profile off
    profile viewer
    keyboard
else
    disp([mfilename,':: exogenous definitions block parsing . ',sprintf('%0.4f',toc),' seconds'])
end

% remove item from block
blocks(current_block_id)=[];
% the equations have been validated, now rebuild them and keep a list of
% the variables defined
DefinedExoList=cell(1,0);%{}
for iii=1:size(ExogenousDefinition_block,1)
    DefinedExoList=[DefinedExoList,ExogenousDefinition_block{iii,1}(1,1)];
    eq_i_='';
    for jjj=1:size(ExogenousDefinition_block{iii,1},2)
        eq_i_=[eq_i_,ExogenousDefinition_block{iii,1}{1,jjj}];
        if ~isempty(ExogenousDefinition_block{iii,1}{2,jjj}) && ExogenousDefinition_block{iii}{2,jjj}~=0
            eq_i_=[eq_i_,'{',sprintf('%0.0f',ExogenousDefinition_block{iii,1}{2,jjj}),'}'];
        end
    end
    ExogenousDefinition_block{iii,1}=eq_i_;
end
% assign information to dictionary
dictionary.exogenous_equations=struct('name',DefinedExoList,'equation',transpose(ExogenousDefinition_block(:,1)));
end
