function h=plot(obj,number,range,varargin)
% PLOT -- plots the transition functions of a VSTAR object
%
% Syntax
% -------
% ::
%
%   h=PLOT(obj,number,range)
%
%   h=PLOT(obj,number,range,varargin)
%
% Inputs
% -------
%
% - **obj** [vstar]: model object
%
% - **number** [integer]: transition function to plot
%
% - **range** [1 x N vector]: discretization of the variable entering the
% transition function
%
% - **varargin** [optional]: pairwise arguments entering matlab's plot
% function
%
% Outputs
% --------
%
% - **h** [handle]:handle to the plot
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:	

g=obj.solution.thresholds{number}.g;

c=num2cell(obj.solution.thresholds{number}.c);

func=obj.thresholds.func;

d=func(range,g,c{:});

h0=plot(range,d,varargin{:});

theName=obj.thresholds(number).name;

if ~strcmp(obj.thresholds(number).lag,'0')
    
    theName=[theName,'{',obj.thresholds(number).lag,'}'];
    
end

title(theName,'interp','none')

if nargout
    
    h=h0;
    
end

end