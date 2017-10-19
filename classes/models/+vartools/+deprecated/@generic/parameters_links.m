function relinkings=parameters_links(obj,estim_names)%,y,s
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


incmnt=300;

% link estimated parameters to the parameter_values matrix
%---------------------------------------------------------
npar=sum(obj.parameters.number);

% links
%------

iter1=0;

nrows1=incmnt;

relinkings=nan(nrows1,2);

nest=numel(estim_names);

for iest=1:nest
    
    est_name=estim_names{iest};
    
    % decompose the name into ploc,chain_loc,state
    %----------------------------------------------
    initialize=iest==1;
    
    [ip,regime_states]=decompose_parameter_name(obj,est_name,initialize);

    % create and store the link between the estimated parameter and the
    % parameter_values matrix
    %---------------------------------------------------------------------
    store_link();

end

relinkings=relinkings(1:iter1,:);

    function store_link()
       
        cols=regime_states;% find(Regimes_mat(:,gov_chain)==istate);
       
        ncols=numel(cols);
       
        if nrows1<iter1+ncols
       
            relinkings(end+incmnt,:)=nan;
       
        end
        
        relinkings(iter1+(1:ncols),:)=[(cols(:)-1)*npar+ip,iest(ones(ncols,1))];
        
        iter1=iter1+ncols;
    
    end

end