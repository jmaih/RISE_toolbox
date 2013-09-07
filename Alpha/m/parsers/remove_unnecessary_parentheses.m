function eqtn=remove_unnecessary_parentheses(eqtn,debug)
% remove parentheses around atoms and multiplicative expressions
% all successive opening and closing parentheses are going down unless they
% are functions. Same for parentheses around multiplicative expressions.
% But in this case, the parentheses cannot be preceeded by / or ^ or
% followed by ^
if nargin<2
    debug=false;
end

char_flag=ischar(eqtn);
if char_flag
    eqtn=cellstr(eqtn);
end

for ieq=1:numel(eqtn)
    eqtn{ieq}=discard_useless_parentheses(eqtn{ieq});
end

if char_flag
    eqtn=char(eqtn);
end

    function eqtn=discard_useless_parentheses(eqtn)
        
        [separators,len,nsep,delimiters]=find_separators(eqtn);
        % legend for the separators matrix
        % 1-type
        % 2-position
        % 3- closing parenthesis, brace or bracket
        % 4- function flag
        % len_delims=length(delimiters);
        discard=false(1,len);
        for isep=1:nsep
            if separators(1,isep)==6 && ~separators(4,isep)
                left=separators(2,isep);
                right=separators(3,isep);
                is_atom=separators(1,isep+1)==7;
                if ~is_atom
                    proceed=left-1==0 || ~any(eqtn(left-1)=='^/');
                    proceed=proceed && (right+1>len || ~any(eqtn(right+1)=='^'));
                    if ~proceed
                        continue
                    end
                    btw_locs= separators(2,:)>left & separators(2,:)<right ;%& separators(1,:)<=len_delims
                    inbetween=delimiters(separators(1,btw_locs));
                    iter=0;
                    depth_parenth=0;
                    depth_brace=0;
                    is_multiplicative=true;
                    while iter<length(inbetween)
                        iter=iter+1;
                        if inbetween(iter)=='('
                            depth_parenth=depth_parenth+1;
                        elseif inbetween(iter)==')'
                            depth_parenth=depth_parenth-1;
                            if depth_parenth<0
                                error(['closing a parenthesis that was not opened in equation ',int2str(ieq)])
                            end
                        elseif inbetween(iter)=='{'
                            depth_brace=depth_brace+1;
                        elseif inbetween(iter)=='}'
                            depth_brace=depth_brace-1;
                            if depth_brace<0
                                error(['closing a parenthesis that was not opened in equation ',int2str(ieq)])
                            end
                        elseif ~depth_parenth && ~depth_brace && any(inbetween(iter)=='+-');
                            if islogical(btw_locs)
                                btw_locs=find(btw_locs);
                            end
                            if separators(2,btw_locs(iter))~=left+1
                                is_multiplicative=false;
                                break
                            end
                        end
                    end
                    if depth_parenth || depth_brace
                        error('open parenthesis not closed')
                    end
                    if debug
                        disp([inbetween,'::',int2str(is_multiplicative)])
                    end
                end
                if is_atom || is_multiplicative
                    discard(separators(2:3,isep))=true;
                end
            end
        end
        %% update equation, separators and discard ...
        eqtn=eqtn(~discard);
    end
end