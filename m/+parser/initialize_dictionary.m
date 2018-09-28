function dictionary = initialize_dictionary()
% INTERNAL FUNCTION
%

dictionary=struct();

dictionary.auxiliary_variables=struct('model',{{}},'ssmodel_solved',{{}});

dictionary.definitions={};

% dictionary.steady_state_parameters={};
dictionary.time_varying_probabilities={};

distr_list=what('+distributions');

for ilist=1:numel(distr_list)

    if ~isempty(distr_list(ilist).m)

        mydistrlist=distr_list(ilist).m;

        break

    end

end
% for some reason I do not understand, this sometimes returns a 2 x 1
% structure instead of a 1 x 1. I have experienced it when using parallel
% computing, but not otherwise.
dictionary.Distributions=strrep(mydistrlist,'.m','_pdf');

dictionary.input_list=parser.input_list();

end


