function Initcond=initial_conditions_to_order_var(Initcond,new_order,options)

% initial_conditions_to_order_var set initial conditions to the order of
% the solution of the model (going from the Alphabetical order)
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


% adjust the start values according to the order_var
%---------------------------------------------------
% Need to take a full otherwise if there are more than 2 dimension,
% indexing does not work.
Initcond.y.y=full(Initcond.y.y); % [variables,ncols,1+n_conditions]
Initcond.y.y=Initcond.y.y(new_order,:,:);% [variables,ncols,1+n_conditions]

iov(new_order)=1:numel(new_order);

% adjust the transition and complementarity functions according to the
% order_var: they will be fed with order_var, but they expect the
% alphabetical order and so the order_var they are fed with has to be
% inverted prior to evaluation.
%--------------------------------------------------------------------------
% Initcond.Qfunc=@(x)Initcond.Qfunc(x(iov));
Initcond.Qfunc=memoize(Initcond.Qfunc,iov);

if ~isempty(Initcond.complementarity)
%     Initcond.complementarity=@(x)Initcond.complementarity(x(iov));
    Initcond.complementarity=memoize(Initcond.complementarity,iov);
end

if ~isempty(Initcond.sep_compl)
%     Initcond.sep_compl=@(x)Initcond.sep_compl(x(iov));
    Initcond.sep_compl=memoize(Initcond.sep_compl,iov);
end

Initcond.simul_honor_constraints_through_switch=options.simul_honor_constraints_through_switch;

Initcond.simul_anticipate_zero=options.simul_anticipate_zero;

end
