function markov_chains_=format_markov_chains(markov_chains_,plist,endo_names)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

nmc=numel(markov_chains_);
cparam=cell(2,nmc+1);
mc_names={markov_chains_.name};
if numel(unique(mc_names))~=numel(mc_names)
    error('two markov chains with the same name')
end
generic_params=regexp(plist,'(a\d+|sig|omg|c)','match');
generic_params=[generic_params{:}];
generic_params=unique(generic_params);
for im=1:nmc
    name=markov_chains_(im).name;
    cparam{1,im}=name;
    controled_params=markov_chains_(im).controled_parameters;
    if isempty(controled_params)
        error(['markov chain "',name,'" must control at least one parameter'])
    end
    % check that the parameters are in the list
    %------------------------------------------
    thisList=[];
    for ilist=1:numel(controled_params)
        thisList=update_sublist(thisList,controled_params{ilist});
    end
    cparam{2,im}=thisList;
    % extract those parameters from the main list
    %--------------------------------------------
    try
        thisListLoc=locate_variables(thisList,plist);
    catch
        disp(thisList)
        error('some parameters in the list above may be controled by different chains')
    end
    plist(thisListLoc)=[];
end
% create the constant-parameter chain and add the remaining
% parameters to the chain
%----------------------------------------------------------
cparam{1,nmc+1}='const';
cparam{2,nmc+1}=plist(:)';

% redo the markov chains
%-----------------------
mc=parser.initialize_markov_chain();

for im=1:size(cparam,2)
    chain_name=cparam{1,im};
    nstates=1;
    pnames=cparam{2,im};
    ptex_names=pnames;
    is_switching=false;
    duration=inf;
    if ~strcmp(chain_name,'const')
        duration=markov_chains_(im).states_expected_duration;
        dr=real(duration);
        di=max(1,imag(duration));
        duration=dr+1i*di;
        nstates=numel(duration);
        if nstates<2
            error(['markov chain "',chain_name,'" should have at least 2 states'])
        end
        is_switching=true;
    end
    mc(im)=parser.initialize_markov_chain(chain_name,nstates,...
        'is_endogenous',false,'param_list',pnames,...
        'param_list_tex',ptex_names,'is_switching',is_switching,...
        'duration',duration);
end
[~,tag]=sort({mc.name});
mc=mc(tag);
markov_chains_=mc(:)'; % for some reason I need this to be a row vector
    function thisList=update_sublist(thisList,cp)
        lp=find(cp=='(');
        middle='';
        if ~isempty(lp)
            rp=find(cp==')');
            middle=cp(lp+1:rp-1);
            cp=cp(1:lp-1);
        end
        if ~any(strcmp(cp,generic_params))
            error(['parameter "',cp,'" unrecognized: maybe you declared a VAR instead of a SVAR?'])
        end
        if ~isempty(middle)
            middle=[middle,','];
            co=find(middle==',');
            prev=0;
            for ico=1:numel(co)
                tmp=middle(prev+1:co(ico)-1);
                if all(isstrprop(tmp,'digit'))
                    % this is the equation
                    eq_loc=str2double(tmp);
                else
                    % this is the equation tag
                    eq_loc=locate_variables(tmp,endo_names);
                end
                something=regexp(plist,[cp,'_',int2str(eq_loc),'_\d+'],'match');
                thisList=union(thisList,[something{:}]);
                prev=co(ico);
            end
        else
            something=regexp(plist,[cp,'\w*'],'match');
            thisList=union(thisList,[something{:}]);
        end
    end
end
