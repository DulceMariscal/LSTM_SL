function update = mySigmoid(e, c1, c2, REV)
    if isempty(REV) || REV==0
        update = (e>=0)./(1+exp(-c1*(e-c2))) - (e<0)./(1+exp(-c1*(-e-c2)));
    else
        update = (e>=0).*(1-1./(1+exp(-c1*(e-c2)))) - (e<0).*(1-1./(1+exp(-c1*(-e-c2))));

%         update = (e>=0).*exp(-c1*(e-c2))./(1+exp(-c1*(e-c2))) - (e<0).*exp(-c1*(-e-c2))./(1+exp(-c1*(-e-c2)));
    end
        
end