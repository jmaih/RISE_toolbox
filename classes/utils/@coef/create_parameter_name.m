function vv=create_parameter_name(vv,endo_names)
% create_parameter_name - creates a parameter name typically for var
% coefficients.
%
% Syntax
% -------
% ::
%
%   vv=create_parameter_name(vv,endo_names)
%
% Inputs
% -------
%
% - **vv** [char|cell]:
%   - char :
%   - cell :
%       - {eqtn,vbl,lag}
%       - {eqtn,vbl,lag,chain,state}
%
% - **endo_names** [cellstr]: list of endogenous variables in the model
%
% Outputs
% --------
%
% - **vv** [char]:
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if ~ischar(vv)% %eqtn,var_pos,lag,chain_name,state
    eqtn=find_variable_position(vv{1});
    var_pos=find_variable_position(vv{2});
    lag=vv{3};
    pname=sprintf('a%0.0f_%0.0f_%0.0f',...
        lag,eqtn,var_pos);
    if numel(vv)>3
        chain_name=vv{4};
        state=vv{5};
        pname=sprintf('%s(%s,%0.0f)',pname,chain_name,state);
    end
    vv=pname;
end
vv=parser.param_texname_to_param_name(vv);
    function loc=find_variable_position(str)
        loc=str;
        if ischar(str)
            if all(isstrprop(str,'digit'))
                warning('expected a variable name but got a string of digit(s)')
                loc=str2double(str);
            else
                loc=find(strcmp(str,endo_names));
                if isempty(loc)
                    error([str,' is not recognized as a variable name'])
                end
            end
        end
        if loc>numel(endo_names)
            disp(vv)
            error(['One of the indexes in the name above ',...
                'exceeds the number of endogenous variables'])
        end
    end
end