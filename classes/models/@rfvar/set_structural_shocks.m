function obj=set_structural_shocks(obj)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 


stshocks=obj.options.structural_shocks;
if ~isempty(stshocks)
    % check that the number of structural shocks does not
    % exceeds the number of exogenous and define the
    % missing structural shocks
    %------------------------------------------------------
    tmp=cell(0,2);
    odd=true;
    for ii=1:numel(stshocks)
        v=stshocks{ii};
        if strcmp(v(1),'"')
            if odd
                error(['two texnames cannot be follow each other and a '...
                    'list cannot start with a texname'])
            end
            tmp{end,2}=v;
        else
            tmp=[tmp;{v,['"',v,'"']}]; %#ok<*AGROW>
        end
        odd=~odd;
    end
    missing=obj.endogenous.number(1)-size(tmp,1);
    if missing<0
        error(['there cannot be more structural shocks than the number of '...
            'reduced-form shocks or endogenous variables'])
    end
    for imiss=1:missing
        v=sprintf('struct_shock_miss_%0.0f',imiss);
        tmp=[tmp;{v,['"',v,'"']}];
    end
    
    % combine with the existing observed exogenous ...
    %---------------------------------------------
    observed_exo_names=get(obj,'exo_list(observed)');
    observed_exo_texnames=strcat('"',get(obj,'exo_tex(observed)'),'"');
    tmp=[tmp
        observed_exo_names(:),observed_exo_texnames(:)];
    stshocks=tmp';
    stshocks=stshocks(:)';
    
    % push the structural shocks
    %---------------------------
    obj=do_names(obj,stshocks,'exogenous');
    
    % redo the observables as some of them could be exogenous
    %--------------------------------------------------------
    obs_names=get(obj,'obs_list');
    obs_texnames=strcat('"',get(obj,'obs_tex'),'"');
    obs=transpose([obs_names(:),obs_texnames(:)]);
    obj=do_names(obj,obs(:)','observables');
end