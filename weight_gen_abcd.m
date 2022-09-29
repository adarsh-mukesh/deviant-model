function [thal_sp, thal_l4_inh, sp_l4_inh, thal_l4_ex, sp_l4_ex, l4_inh_ex]=weight_gen_abcd(sigma)

global sp n1 n2 n3
weight_init=0.2*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*sigma^2));
for i=1:sp
    thal_sp(i,i:i+2*n1)=weight_init;
end
%-------------------------------------
l4_inh=sp-2*n2;
l4_ex=l4_inh-2*n3;
%----------------------------------------------------
%----------------------------------------
%connections from thalamic to layer 4 (flat)
for i=1:l4_inh
    for j=1:2*n1+1
       thal_l4_inh(i,j+i-1)=0.02;
    end
    sp_l4_inh(i,i:i+2*n2)=0.11;
    %thal_l4(i,i:i+8)=2*weight_init;
end
%-------------------------------------------
for i=1:l4_ex
    for j=1:2*n1+1
       thal_l4_ex(i,j+i-1)=0.02;
    end
    sp_l4_ex(i,i:i+2*n2)=0.11;
    %thal_l4(i,i:i+8)=2*weight_init;
    l4_inh_ex(i,i:i+2*n3)=0.06;
end
