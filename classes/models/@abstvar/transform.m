function p=transform(self,p)

linres=self.estim_.linres;

if ~isempty(linres)
    
    p=linres.a_to_a2tilde(p);
    
end

end