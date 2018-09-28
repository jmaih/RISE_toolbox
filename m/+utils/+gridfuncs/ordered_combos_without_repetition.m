function Sols=ordered_combos_without_repetition(choice_set,m)

% choice_set: set to choose from
% m: number of elements to choose

n=numel(choice_set);

if m>n
    
    error('m cannot be greater than the # setx')
    
end

protopath=nan(1,m);

siz_sols=utils.gridfuncs.my_nchoosek(n,m);

Sols=zeros(siz_sols,m);

itersol=0;

follow_path(protopath,1,1:n)

if itersol~=siz_sols
    
    warning('unpleasant arithmetics')
    
    Sols=Sols(1:itersol,:);
    
end

Sols=choice_set(Sols);

    function follow_path(path1,next_pos,xchoice_set)
        
        if nargin<3
            
            next_pos=1;
            
        end
        
        while ~isempty(xchoice_set) % && (nchoice>=nancount)
            
            this_guy=xchoice_set(1);
            
            path1(next_pos)=this_guy;
            
            if next_pos==m
                
                add_sol(path1)
                
            end
            
            xchoice_set=xchoice_set(2:end);
            
            follow_path(path1,next_pos+1,xchoice_set);
            
        end
        
    end

    function add_sol(path1)
        
        itersol=itersol+1;
        
        Sols(itersol,:)=path1;
        
    end

end

%{ 
ALTERNATIVE ALGORITHM
function Sols=all_combos_without_repetition(setx,m)

n=numel(setx);

if m>n
    
    error('m cannot be greater than the # setx')
    
end

protopath=nan(1,m);

incr=1000;

siz_sols=incr;

Sols=zeros(siz_sols,m);

itersol=0;

follow_path(1,protopath)

Sols=Sols(1:itersol,:);

Sols=setx(Sols);

    function follow_path(start,path1,next_pos)
        
        if nargin<3
            
            next_pos=1;
            
        end
        
        for ii=start:n
            
            next_start=ii+1;
            
            % accept
            %--------
            path1(next_pos)=ii;
            
            if next_pos==m
                
                add_sol(path1)
                
            end
            
            follow_path(next_start,path1,next_pos+1);
            
        end
        
    end

    function add_sol(path1)
        
        itersol=itersol+1;
        
        if itersol>siz_sols
            
            Sols=[Sols
                zeros(incr,m)];
            
            siz_sols=siz_sols+incr;
            
        end
        
        Sols(itersol,:)=path1;
        
    end

end
%}