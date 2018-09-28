function strcell=poly2str(Pset,varnames)
% POLY2STR - convert multivariate polynomials from arrays to strings
% usage: strcell=poly2str(Pset,varnames)
%
% INPUTS: Pset, varnames
%  Pset is a cell array of polynomial coefficients where row j of Pset{i}
%    is [c,k1,k2,...] representing the monomial c*x1^k1*x2^k2*...
%    Pset can be the output of petschel.groebner(...)
%  varnames (optional) is a cell array of strings of length at least the
%    maximum dimension of the arrays in Pset [default: {'x1','x2',...}]
%
% OUTPUTS: strcell
%  strcell is a cell array of strings, one per polynomial
%
% EXAMPLE:
%  poly2str({[1,1,2]}) % returns {'x1*x2^2'}
%
% SEE ALSO:
%  petschel.groebner, petschel.str2poly

% Author: Ben Petschel 19/6/2009
%
% Change history:
%  19/6/2009 - first release
%  20/3/2010 - changed array representation to rectangular instead of n-d
%  2/11/2010 - bugfix for maxnvars calc that was broken in the array rep update

maxnvars = max(cellfun(@(x)size(x,2)-1,Pset)); % determine nvars from array size
if nargin<2,
  varnames = cell(1,maxnvars);
  for i=1:maxnvars,
    varnames{i} = sprintf('x%d',i);
  end;
elseif maxnvars>numel(varnames),
  error('not enough variable names provided');
end;

strcell = cell(size(Pset));

for i=1:numel(Pset),
  str = '';
  P = Pset{i};
  s = size(P);
  c = cell(1,length(s));
  for j=size(P,1):-1:1,
    if P(j,1)~=0,
      if (P(j,1)>0) && ~isempty(str),
        str = [str,'+'];
      end;
      if P(j,1)==1,
        last1 = true;
      elseif P(j,1)==-1,
        % only print the '-'
        last1 = true;
        str = [str,'-'];
      else
        last1 = false;
        str = [str,sprintf('%g',P(j,1))]; % include coefficient
      end;
      deg = P(j,2:end);
      for k=1:length(deg),
        if deg(k)>0,
          if last1,
            last1 = false;
          else
            str = [str,'*'];
          end;
          str = [str,varnames{k}]; % print variable name
        end;
        if deg(k)>1,
          str = [str,sprintf('^%d',deg(k))]; % print exponent
        end;
      end;
      if last1,
        % this must be a constant term 1, so print the 1
        str = [str,'1'];
      end;
    end;
  end;
  strcell{i}=str;
end;

end % main function poly2str(...)
