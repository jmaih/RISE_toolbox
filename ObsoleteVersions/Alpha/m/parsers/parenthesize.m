function eqtn_out=parenthesize(eqtn_in,func,debug)
% func is any of the following: '^', '/' , '-', '+' ,'*'
% this function is called by sad_forward, rise_algebra_cat

if nargin<3
    debug=false;
end

char_flag=ischar(eqtn_in);
if char_flag
    eqtn_in=cellstr(eqtn_in);
end
eqtn_out=eqtn_in;
for ii=1:numel(eqtn_in)
    eqtn_out{ii}=parenthesize_intern(eqtn_in{ii});
end
if char_flag
    eqtn_out=char(eqtn_out);
end

    function eqtn=parenthesize_intern(eqtn)
        length_string=length(eqtn);
        
        mult_ready=true;
        power_ready=true;
        div_ready=true;
        minus_ready=true;
        plus_ready=true;
        ready = (func=='*' && mult_ready)||...
            (func=='^' && power_ready)||...
            (func=='/' && div_ready)||...
            (func=='-' && minus_ready)||...
            (func=='+' && plus_ready);
        if ~strcmp(func,'+')
            iter=0;
            while iter<length_string  && ready
                iter=iter+1;
                token=eqtn(iter);
                if token=='('
                    % try and close the parenthesis before proceeding
                    depth=1;
                    while depth && iter<length_string
                        iter=iter+1;
                        token=eqtn(iter);
                        if token=='('
                            depth=depth+1;
                        elseif token==')'
                            depth=depth-1;
                        end
                    end
                    if depth
                        error('parentheses not closed')
                    end
                elseif  token=='+'
                    mult_ready=false;
                    power_ready=false;
                    div_ready=false;
                    minus_ready=false;
                elseif  token=='-'
                    mult_ready=false;
                    power_ready=false;
                    div_ready=false;
                    minus_ready=false;
                elseif  token=='/'
                    power_ready=false;
                elseif  token=='*'
                    power_ready=false;
                elseif  token=='^'
                    power_ready=false;
                elseif token==')'
                    error(['parenthesis at position ',int2str(iter),' closing before an opening'])
                end
                ready = (func=='*' && mult_ready)||(func=='^' && power_ready)||...
                    (func=='/' && div_ready)||(func=='-' && minus_ready);
                % Note that at this stage I do not check the plus since it is
                % always ready to go.
            end
        end
        
        if debug
            fprintf(1,'mult_ready=%0.1g, power_ready=%0.1g, div_ready=%0.1g, minus_ready=%0.1g , ratio=%0.3g percent \n',...
                mult_ready,power_ready,div_ready,minus_ready,100*iter/length_string);
        end
        
        if ~ready
            eqtn=strcat('(',eqtn,')');
        end
    end
end

