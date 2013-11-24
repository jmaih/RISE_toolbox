function wrt=extend_differentiation_list(incidence,wrt0)

set_differentiation_list();

incidence=incidence(:);
n=numel(incidence);
wrt=cell(1,n);
wrt(incidence)=wrt0;
iter=0;
for ii=1:n
    if isempty(wrt{ii})
        iter=iter+1;
        wrt{ii}=['zZzZzZz_',sprintf('%0.0f',iter)];
    end
end

    function set_differentiation_list()
        reprocess=size(wrt0,2)==2 && isnumeric(wrt0{1,2});
        if reprocess
            with_respect_to=cell(1,300);
            jter=0;
            for iii=1:size(wrt0,1)
                digits=wrt0{iii,2};
                xx=wrt0{iii,1};
                for id=1:numel(digits)
                    jter=jter+1;
                    if jter==size(with_respect_to,2)
                        with_respect_to{end+300}={}; %#ok<AGROW>
                    end
                    with_respect_to{jter}=[xx,'_',sprintf('%0.0f',digits(id))]; 
                end
            end
            wrt0=with_respect_to(1:jter);
        end
    end
end