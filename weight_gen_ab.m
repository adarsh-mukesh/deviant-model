function [thal_sp]=weight_gen_ab(sigma)
global sp n1
weight_init=0.1*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*sigma^2));
for i=1:sp
    thal_sp(i,i:i+2*n1)=weight_init;
end
end