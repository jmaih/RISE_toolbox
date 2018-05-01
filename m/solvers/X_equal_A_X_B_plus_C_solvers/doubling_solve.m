function [P,retcode,good]=doubling_solve(A,B,C,options)
% doubling_solve solves the linear equation X=A*X*B+C
%
% ::
%
%   [P,retcode]=doubling_solve(A,B,C)
%   [P,retcode]=doubling_solve(A,B,C,options)
%
% Args:
%    - A :
%    - B :
%    - C :
%    - options :
%
% Returns:
%    :
%    - P :
%    - retcode :
%
% Note:
%
% Example:
%
%    See also:

mydefaults=the_defaults();

if nargin==0
    
    if nargout
        
        P=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

if nargin<4
    
    options=[];
    
    if nargin<3
        
        error([mfilename,':: at least 3 arguments should be provided'])
    
    end
    
end

default_solve=disp_defaults(mydefaults);

if isempty(B),B=A';end

if isempty(A),A=B';end

if isempty(options)
    
    options=default_solve;
    
else
    
    fnames=fieldnames(default_solve);
    
    for itype=1:size(mydefaults,1)
        
        name=fnames{itype};
        
        if ~isfield(options,name)
            
            options.(name)=default_solve.(name);
            
        end
        
    end
    
end

P0=C;

symmetric=isequal(A,B');

Gl=A;

if ~symmetric
    
    Gr=B;
    
end

[P,~,retcode]=fix_point_iterator(@iterator,P0,options);

good=true(size(P,1),1);

    function [P,F0]=iterator(P0)
        
        if symmetric
            
            P=P0+Gl*P0*Gl';
            
        else
            
            P=P0+Gl*P0*Gr;
            
            Gr=Gr*Gr;
            
        end
        
        Gl=Gl*Gl;
        
        F0=P-P0;
        
    end

end

function d=the_defaults()

d=cell(0,4);

end