function eqtn_out=remove_zeros_and_ones(eqtn_in,debug)
if nargin<2
    debug=false;
end
cellmode=iscellstr(eqtn_in);
if ~cellmode
    eqtn_in=cellstr(eqtn_in);
end
eqtn_out=eqtn_in;
neq=numel(eqtn_out);
if neq>1
    debug=false;
end
%% zeros in loops
for ieq=1:neq
    eqtn_out{ieq}=remove_zeros_and_ones_intern(eqtn_out{ieq});
end

%% redo if necessary
bad=~cellfun(@isequal,eqtn_in,eqtn_out);
if any(any(bad))
    eqtn_out(bad)=string_optimize.remove_zeros_and_ones(eqtn_out(bad),debug);
end
if ~cellmode
    eqtn_out=char(eqtn_out);
end

%% nested functions
    function eqtn_out=remove_zeros_and_ones_intern(eqtn_in)
        eqtn_out=eqtn_in;
        if debug
            disp('original equation')
            disp(eqtn_out)
        end
        locs=regexp(eqtn_out,'(?<![\w.])0(?![.\w\^])','start'); % <-- locs=find(eqtn=='0');
        if isempty(locs)
            return
        end
        eqtn_length=length(eqtn_out);
        discard=false(1,eqtn_length);
        for ii=1:numel(locs)
            if discard(locs(ii))
                continue
            end
            % look around
            left='';
            right='';
            if locs(ii)>1
                left=eqtn_out(locs(ii)-1);
            end
            if locs(ii)<eqtn_length
                right=eqtn_out(locs(ii)+1);
            end
            
            if any(strcmp(left,{'/',')'}))
                error('')
            end
            if any(strcmp(right,{'('}))
                error('')
            end
            if strcmp(left,'-')
                left='+';
            end
            switch [left,'0',right]
                case '0'
                    return
                case {'0)','(0'}
                    error('')
                case {'+0','+0+','+0-','+0)'}
                    % discard both. If the first case applies and the end
                    % string is empty, then it will just be replace by 0
                    discard([locs(ii)-1,locs(ii)])=true;
                case {'0+','0-'}
                    discard(locs(ii))=true;
                case {'+0*','+0/','(0*','(0/','0*','0/'}
                    start=locs(ii)-strcmp(left,'+');
                    iter=parse_left_right(start,'r');
                    discard(start:iter)=true;
                    % take action to the right and then see if you can discard
                    % further
                case {'^0','^0+','^0-','^0)'}
                    % take action to the left and see if you can discard
                    start=locs(ii);
                    iter=parse_left_right(start,'l');
                    discard(iter:start-1)=true;
                    eqtn_out(locs(ii))='1';
                case {'*0*','*0/','^0*','^0/'}
                    start=locs(ii);
                    iter_l=parse_left_right(start,'l');
                    if strcmp(left,'^')
                        eqtn_out(locs(ii))='1';
                        discard(iter_l:start-1)=true;
                    else
                        iter_r=parse_left_right(start,'r');
                        discard(iter_l:iter_r)=true;
                    end
                    % take action to the right and to the left
                case {'*0','*0+','*0-','*0)'}
                    start=locs(ii);
                    iter=parse_left_right(start,'l');
                    discard(iter:start)=true;
            end
            if debug
                disp(['case ',left,'0',right])
                disp(eqtn_out(~discard))
            end
        end
        
        %% Remove any "1^": if not preceeded by a word character or a .
        locs1=regexp(eqtn_out,'(?<![\w\.])1\^','start');
        for ii=1:numel(locs1)
            if discard(locs1(ii))
                continue
            end
            start=locs1(ii)+1;
            iter=parse_left_right(start,'r');
            discard(start:iter)=true;
        end
        %% now discarding
        eqtn_out(discard)=[];
        eqtn_out=strrep(eqtn_out,'()','0');
        if isempty(eqtn_out)
            eqtn_out='0';
        end
        if ~strcmp(eqtn_in,eqtn_out)
            eqtn_out=remove_zeros_and_ones_intern(eqtn_out);
        end
        
        function iter=parse_left_right(start,lr)
            if strcmp(lr,'r')
                left_='(';
                right_=')';
                incmnt=1;
                checkfun=@(a)lt(a,eqtn_length);
            elseif strcmp(lr,'l')
                left_=')';
                right_='(';
                incmnt=-1;
                checkfun=@(a)gt(a,2);
            end
            iter=start;
            depth=0;
            while checkfun(iter)
                iter=iter+incmnt;
                if strcmp(eqtn_out(iter),left_)
                    depth=depth+1;
                elseif strcmp(eqtn_out(iter),right_)
                    depth=depth-1;
                    % the parenthesis was probably opened long before/after the
                    % start. so exit
                    if depth<0;
                        % we've gone too far. correct that real quick
                        iter=iter-incmnt;
                        break
                    end
                elseif depth==0 && any(eqtn_out(iter)=='+-')
                    % we've gone too far. correct that real quick
                    iter=iter-incmnt;
                    break
                else
                    % if we come to the end of the string, there is no
                    % correction needed
                end
            end
        end
    end
end