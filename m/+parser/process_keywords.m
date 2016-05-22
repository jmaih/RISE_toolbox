function xin=process_keywords(xin,endo_names)

iii_=[]; 

n=[];

myreplace=@replacer; %#ok<NASGU>

kwords={'movave','movavg','movgeom','movprod','movsum','diff','difflog','dot'};

kwords=parser.cell2matize(kwords);

trigger_patt=['\<',kwords,'\>\('];

endo_names_=parser.cell2matize(endo_names);

newpatt_0=['\<',endo_names_,'\>'];

newpatt_d=['\<',endo_names_,'\>(\(|\{)((\+|\-)?\d+)(\)|\})'];

wedge='_____';

for i1=1:size(xin,1)
    
    rawline=xin{i1,2};
    
    file_name=xin{i1,3};
    
    line_number=xin{i1,1};
    
    while 1
        
        [start_index,end_index]=regexp(rawline,trigger_patt,'start','end');
        
        if isempty(start_index)
            
            break
            
        end
        % process last only
        start_index=start_index(end);
        
        end_index=end_index(end); 
        
        depth=0;
        
        left=rawline(1:start_index-1);
        
        right=rawline(end_index+1:end);
        
        start=find_closing();
        
        stud=rawline(start_index:end_index-1);
        
        middle=right(1:start-1);
        
        right=right(start+1:end);
        
        comma=find(middle==',',1,'last');
        
        main_matter=middle;
        
        process_matter();
        
        rawline=[left,main_matter,right];
        
        xin{i1,2}=rawline;
        
    end
    
end

    function c=find_closing()
        
        closingFound=false;
        
        c=0;
        
        while ~closingFound && c<length(right)
            
            c=c+1;
            
            if right(c)=='('
                
                depth=depth+1;
                
            elseif right(c)==')'
                
                closingFound=depth==0;
                
                if ~closingFound
                    
                    depth=depth-1;
                    
                end
                
            end
            
        end
        
        if ~closingFound
            
            error(['Missing closing for pseudo function in ',...
                 file_name,' at line ',sprintf('%0.0f',line_number)])
            
        end
        
    end

    function process_matter()
        
        n=[];
        
        if ~isempty(comma)
            
            nstar=eval(middle(comma+1:end));
            
            if ~isnan(nstar)
                
                main_matter=middle(1:comma-1);
                
                n=nstar;
                
            end
            
            if ~isfinite(n)||(floor(n)~=ceil(n))
                
                error(['wrong specification of integer when using ',stud,...
                     'in ',file_name,' at line ',sprintf('%0.0f',line_number)])
                
            end
            
        end
        
        v=['(',main_matter,')'];
        
        switch stud
            
            case {'movave','movavg'}
                
                if isempty(n),n=-4;end
                
                str=many_together('+');
                
                str=['1/',int2str(abs(n)),'*',str];
                
            case 'movprod'
                
                if isempty(n),n=-4;end
                
                str=many_together('*');
                
            case 'movgeom'
                
                if isempty(n),n=-4;end
                
                str=many_together('*');
                
                str=[str,'^(1/',int2str(abs(n)),')'];
                
            case 'movsum'
                
                if isempty(n),n=-4;end
                
                str=many_together('+');
                
            case 'diff'
                
                if isempty(n),n=-1;end
                
                str=[v,'-',apply_lag(v,abs(n))];
                
            case 'difflog'
                
                if isempty(n),n=-1;end
                
                str=['log(',v,')-log(',apply_lag(v,abs(n)),')'];
                
            case 'dot'
                
                if isempty(n),n=-1;end
                
                str=[v,'/',apply_lag(v,abs(n))];
                
            otherwise
                
                error(['Problem in ',file_name,' at line ',...
                    sprintf('%0.0f',line_number),' Please report this ',...
                    'error to the forum or contact junior.maih@gmail.com'])
                
        end
        
        main_matter=['(',str,')'];
        
        function str=many_together(sign_)
            
            str=['(',v];
            
            for ii=1:abs(n)-1
                
                str=[str,sign_,apply_lag(v,ii)]; %#ok<AGROW>
                
            end
            
            str=[str,')'];
            
        end
        
    end

    function v=apply_lag(v,ii)
        
        iii_=ii*sign(n);
        % Apply leads and lags: process 0 and +-n separately. Process the
        % guys with lead or lags first!
        %-----------------------------------------------------
        v=regexprep(v,newpatt_d,'${myreplace($1,$3)}');
        
        v=regexprep(v,newpatt_0,'${myreplace($1)}');
        
        v=strrep(v,wedge,'');
        
        v=strrep(v,'{0}','');
        
    end
        
    function str=replacer(v,d)
        
        if nargin<2
            
            d='0';
            
        end
        
        d=str2double(d);
        
        if ~isfinite(d)||(floor(d)~=ceil(d))
            
            error(['wrong specification of integer when using ',stud,...
                 'in ',file_name,' at line ',sprintf('%0.0f',line_number)])
            
        end
        
        d=int2str(d+iii_);
        
        str=[v,wedge,'{',d,'}'];
        
    end
    
end

%--------------------------------------------------------------------------
