function Pset=str2poly(strcell,varnames)
% POLY2STR - convert multivariate polynomials from strings to arrays
% usage: Pset=petschel.poly2str(strcell,varnames)
%
% INPUTS: strcell, varnames
%  strcell is a cell array of strings, one per polynomial
%    all numbers must be whole or decimal form (no scientific notation)
%    brackets and double negatives are not allowed.
%    All spaces in the strings are ignored.
%  varnames (optional) is a cell array of strings which must include all
%    variable names that occur in strcell [default: {'x1','x2',...}]
%
% OUTPUTS: Pset
%  Pset is a cell array of polynomial coefficients where row j of Pset{i}
%    is [c,k1,k2,...] representing the monomial c*x1^k1*x2^k2*...
%
% EXAMPLE:
%  str2poly({'x1*x2^2'})    % returns {[1,1,2]}
%  str2poly({'-x1*x2^2'})   % returns {[-1,1,2]}
%  str2poly({'-1*x1*x2^2'}) % returns {[-1,1,2]}
%  str2poly({'- 1 * x 1 * x 2 ^ 2'}) % returns {[-1,1,2]}
%
% SEE ALSO:
%  petschel.groebner, petschel.poly2str

% Author: Ben Petschel 19/6/2009
%
% Change history:
%  19/6/2009 - first release
%  22/6/2009 - fixed bug in handling "-x1" etc.
%  17/7/2009 - allow whitespace in strings
%  20/3/2010 - changed array representation to rectangular instead of n-d

if nargin<2
    
    varnames = {}; % will fill in varnames later, as needed
    
end

Pset = cell(size(strcell));

for i=1:numel(strcell)
    
    str = strcell{i};
    
    % ignore all spaces
    remain = str;
    
    remain = strrep(remain,' ','');
    
    % first insert "+" before every "-", for easier term separation later
    str = '';
    
    while ~isempty(remain)
        
        if remain(1)=='-'
            
            if isempty(str)
                
                str = '-';
                
            else
                
                str = [str,'+-']; %#ok<*AGROW>
                
            end
            
        end
        
        [tok,remain]=strtok(remain,'-'); %#ok<*STTOK>
        
        str = [str,tok];
        
    end
    
    % now separate into terms
    remain = str;
    
    terms = [];
    
    while ~isempty(remain)
        
        coeff = 1;
        
        exps = [];
        
        [term,remain]=strtok(remain,'+');
        
        % process term, separating into products
        rem2 = term;
        
        while ~isempty(rem2)
            
            [part,rem2]=strtok(rem2,'*');
            
            if ~isnan(str2double(part))
                % part is a number, multiply it by the coefficient
                coeff = coeff * str2double(part);
                
            else
                % part is var^n or var or -var^n or -var
                if part(1)=='-'
                    
                    coeff = -coeff;
                    
                    part=part(2:end);
                    
                    if isempty(part)
                        
                        error('isolated minus sign');
                        
                    end
                    
                end
                
                ind = find(part=='^');
                
                if isempty(ind)
                    
                    var = part;
                    
                    expo = 1;
                    
                else
                    
                    ind = ind(1);
                    
                    var = part(1:ind-1);
                    
                    expo = str2double(part(ind+1:end));
                    
                    if isnan(expo)
                        % exponent is not valid
                        error('invalid exponent in part %s',part);
                        
                    end
                    
                end
                
                if isempty(varnames)
                    % see if var is of form 'x1' etc
                    if (length(var)>1) && (var(1)=='x') && ~isnan(str2double(var(2:end)))
                        
                        % variable name is ok
                        k = str2double(var(2:end));
                        
                    else
                        
                        error('variable name not of form "x1","x2",...');
                        
                    end
                    
                else
                    
                    ind = find(strcmp(var,varnames));
                    
                    if isempty(ind)
                        
                        error('variable name not in list');
                        
                    end
                    
                    k = ind(1);
                    
                end
                
                if length(exps)<k
                    
                    exps = [exps,zeros(1,k-length(exps))];
                    
                end
                
                exps(k) = exps(k)+expo;
                
            end
            
        end % while ~isempty(rem2),
        
        % term is processed, now add it to the list of terms
        if size(terms,2)<(length(exps)+1)
            
            terms = [terms,zeros(size(terms,1),length(exps)+1-size(terms,2))];
            
        elseif length(exps)<size(terms,2)-1
            
            exps = [exps,zeros(1,size(terms,2)-1-length(exps))];
            
        end
        
        terms = [terms;coeff,exps];
        
    end % while ~isempty(remain),
    
    Pset{i}=terms;
    
end % for i=1:numel(strcell),


end % main function str2poly(...)
