function str=optimize(str)

% replaces

str_length0 = length(str);
str_length = str_length0;

parentheses=[];
get_parentheses();

optimize_sign();

myset=[1,2,4];%,3,5

if any(myset==1)
    % replace '((....))'  with '(....)'
    lplp_loc=strfind(str,'((');
    % remove the overlaps
    for ii=2:numel(lplp_loc)
        if lplp_loc(ii)==lplp_loc(ii-1)+1
            lplp_loc(ii)=nan;
        end
    end
    lplp_loc(isnan(lplp_loc))=[];
    dejavu=[];
    for iloc=1:numel(lplp_loc)
        nail_it=find(parentheses(1,:)==lplp_loc(iloc)); % take one at a time
        ll=parentheses(1,nail_it);
        rr=parentheses(2,nail_it);
        if strcmp(str(rr-1),')') && (ll==1||any(str(ll-1)=='+-*/'))
            nail_close=find(parentheses(2,:)==rr-1);
            if nail_close==nail_it+1 % penultimate closing =? post nail
                dejavu=[dejavu,ll,rr]; %#ok<*AGROW>
            end
        end
    end
    delete_string(dejavu);
end

% replace (n) n
optimize_digits_in_parentheses();

if any(myset==2)
    optimize_1();
    % remove any 1^ --->1
    dejavu=[];
    powerOne_loc=strfind(str,'1^');
    for iloc=1:numel(powerOne_loc)
        ll=powerOne_loc(iloc)+1;
        if strcmp(str(ll+1),'(')
            par_loc= parentheses(1,:)==ll+1;
            rr=parentheses(2,par_loc);
            dejavu=[dejavu,(ll:rr)];
        else% digits...
            rr=ll+1;
            while rr+1<str_length && any(str(rr+1)=='.0123456789')
                rr=rr+1;
            end
            dejavu=[dejavu,(ll:rr)];
        end
    end
    delete_string(unique(dejavu));
end


if any(myset==3)
    optimize_0();
    % ()^0 ---> 1 n^0  ---> 1
    dejavu=[];
    pow_zero=strfind(str,'^0');
    for iz=1:numel(pow_zero)
        rr=pow_zero(iz);
        if rr+1==str_length || any(str(rr+2)=='+-/*')
            str(rr+1)='1';
            % now search for ll backwards
            ll=rr-1;
            if strcmp(str(ll),')')
                par_loc=parentheses(2,:)==ll;
                opening=parentheses(1,par_loc);
                if opening==1||~isletter(str(opening-1))
                    ll=opening;
                    dejavu=[dejavu,(ll:rr)];
                end
            else % digits ? could well be a function without arguments...
                while ll-1>1 && any(str(ll-1)=='.0123456789')
                    ll=ll-1;
                end
                dejavu=[dejavu,(ll:rr)];
            end
        end
    end
    delete_string(dejavu);
end

if any(myset==4)
    % remove unneeded parentheses
    optimize_parentheses_around_several_atoms();
% % % % %     % (atom)--->atom
% % % % %     % {'',+}(atom+atom+...+atom){-+,''}--->atom
% % % % %     dejavu=[];
% % % % %     for ipar=1:size(parentheses,2)
% % % % %         ll=parentheses(1,ipar);
% % % % %         rr=parentheses(2,ipar);
% % % % %         if rr>ll+1
% % % % %             if (is_atom(str(ll+1:rr-1)) && (ll==1||isequal(str(ll-1),'+')))||...
% % % % %                     ((ll==1 ||isequal(str(ll-1),'+')) && ... %((ll==1 ||~ismember(str(ll-1),'^,*-/')) && ...
% % % % %                     (rr==str_length ||~any(str(rr+1)=='^*/')))
% % % % %                 dejavu=[dejavu,ll,rr];
% % % % %             end
% % % % %         end
% % % % %     end
% % % % %     delete_string(dejavu);
end

if any(myset==5)
    
    % replace 0*() and 0/() with 0
    dejavu=[];
    zero_times=strfind(str,'0*');
    zero_times=union(zero_times,strfind(str,'0/'));
    for iz=1:numel(zero_times)
        ll=zero_times(iz)+1;
        if ~ismember(ll,dejavu)
            if ll==2||any(str(ll-2)=='(+-*')
                rr=ll+1;
                if strcmp(str(rr),'(')
                    nail=parentheses(1,:)==rr;
                    rr=parentheses(2,nail);
                elseif any(str(rr)=='.0123456789')
                    while rr+1<str_length && any(str(rr+1)=='.0123456789')
                        rr=rr+1;
                    end
                else
                    % perhaps it is a variable
                    [~,rest]=strtok(str(rr:end),'+-*/^()');
                    if isempty(rest)
                        rr=str_length;
                    elseif any(rest(1)=='+-*/^)') % a variable
                        while rr+1<str_length && ~any(str(rr+1)=='+-*/^)')
                            rr=rr+1;
                        end
                    else % a function
                        nail_par=find(parentheses(1,:)>rr,1,'first');
                        if ~isempty(nail_par)
                            rr=parentheses(2,nail_par);
                        end
                    end
                    % it is a variable or a function, look for the next opening
                    % parenthesis
                end
                dejavu=[dejavu,(ll:rr)];
            end
        end
    end
    delete_string(unique(dejavu));
    
    % replace 0^() with 0 ... dangerous if the exponent is 0
    dejavu=[];
    zero_pow=strfind(str,'0^');
    for iz=1:numel(zero_pow)
        ll=zero_pow(iz)+1;
        rr=ll+1;
        if any(str(rr)=='.0123456789')
            while rr+1<str_length && any(str(rr+1)=='.0123456789')
                rr=rr+1;
            end
            if strcmp(str(ll+1:rr),'0')
                str(ll-1)='1';
            end
            dejavu=[dejavu,(ll:rr)];
        end
    end
    delete_string(unique(dejavu));
end

if str_length ~= str_length0
    str=rise_sad.optimize(str);
end

    function delete_string(indexes)
        if ~isempty(indexes)
            str(indexes)=[];
            str_length=length(str);
            % there should be a more efficient way of updating parentheses...
            get_parentheses();
        end
    end

    function get_parentheses()
        right=strfind(str,')');
        left=strfind(str,'(');
        parentheses=[left;right];
        for ip=numel(left):-1:1
            closing=find(right>left(ip),1,'first');
            parentheses(2,ip)=right(closing);
            right(closing)=[];
        end
    end

    function optimize_digits_in_parentheses()
        str=regexprep(str,'(?<![\w])(\()(\d*\.?\d*)(\))','$2'); % <--- str=regexprep(str,'(\W)(\()(\d*\.?\d*)(\))','$1$3');
        str=regexprep(str,'(?<![\w])(\()(y|p|param|x|def|ss)(\_)(\d*)(\))','$2$3$4');
        get_parentheses();
        str_length=length(str);
    end

    function optimize_sign()
        str=strrep(str,'+-','-');
        str=strrep(str,'-+','-');
        str=strrep(str,'++','+');
        str=strrep(str,'--','+');
        get_parentheses();
        str_length=length(str);
    end

    function optimize_1()
        % 1* cannot be preceded by any [a_zA_Z0_9./]
        str=regexprep(str,'(?<![\w|/|.])1\*',''); 
        % *1 cannot be followed by any [a_zA_Z0_9./]
        str=regexprep(str,'\*1(?![\w|/|.])',''); 
        % ^1 cannot be followed by any [a_zA_Z0_9./]
        str=regexprep(str,'\^1(?![\w|/|.])',''); 
        get_parentheses();
        str_length=length(str);
    end

    function optimize_0()
        % 0- or 0+ cannot be preceded by any [a_zA_Z0_9./]
        str=regexprep(str,'(?<![\w|/|.])(0[\-|\+])',''); 
        % -0 cannot be followed by any [a_zA_Z0_9./]
        str=regexprep(str,'(\-0)(?![\w|/|.])',''); 
        % +0 cannot be followed by any [a_zA_Z0_9./]
        str=regexprep(str,'(\+0)(?![\w|/|.])',''); 
        % () cannot be preceded by any [a_zA_Z0_9]: These may have been created by the removal of -0 +0 0+ 
        str=regexprep(str,'(?<![\w])(\(\))',''); 
        get_parentheses();
        str_length=length(str);
    end
    function optimize_parentheses_around_several_atoms()
        % (atom+atom*atom/atom^atom) cannot be preceded by any [\w^*-] and
        % cannot be followed by [*/]
        NotBefore='(?<![-*/^\w])';
        leftPar='(\()';
        rightPar='(\))';
        Middle='([\w*+^()]*)';
        NotAfter='(?![*/])';
        pat=[NotBefore,leftPar,Middle,rightPar,NotAfter];
        str=regexprep(str,pat,'$2'); 
        get_parentheses();
        str_length=length(str);
    end

end
