function [retcode,G0,G1,G2]=dynamic_derivatives(deriv_code,y,x,ss,param,def,order)

%anonymous function for checking validity
%----------------------------------------
valid=@(x)~any(isnan(x(:))) && ~any(isinf(x(:))); % nans in jacobian

h=size(param,2);
hy=size(y,2); % if the steady state is unique, hy may be 1 while h >1
if ~ismember(hy,[1,h])||size(x,2)~=hy||size(ss,2)~=hy
    error('columns of y,x,ss inconsistent with the parameters')
end

% initialize output
%------------------
G0=cell(h);
G1=cell(h);
G2=cell(h);

def01=[];
for s0=1:h
    for s1=1:h
        y01=y(:,min(hy,s1));
        x01=x(:,min(hy,s1));
        ss01=ss(:,min(hy,s1));
        if ~isempty(def)
            def01=def(:,min(hy,s1));
        end
        [retcode,G0{s0,s1},G1{s0,s1},G2{s0,s1}]=evaluate_derivatives(s0,s1);
        if retcode
            return
        end
    end
end

    function [retcode,G0,G1,G2]=evaluate_derivatives(s0,s1)
        G0=[];G1=[];G2=[];
       switch order
            case 0
                G0=online_function_evaluator(deriv_code,y01,x01,ss01,param(:,s1),def01,s0,s1);
            case 1
                [G0,G1]=online_function_evaluator(deriv_code,y01,x01,ss01,param(:,s1),def01,s0,s1);
            case 2
                [G0,G1,G2]=online_function_evaluator(deriv_code,y01,x01,ss01,param(:,s1),def01,s0,s1);
            otherwise
       end
        retcode=2*(~valid(G0)||(order>0 && ~valid(G1))||(order>1 && ~valid(G2)));
    end
end