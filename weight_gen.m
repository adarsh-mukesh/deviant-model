function [thal_sp, thal_l4, sp_l4] = weight_gen(sp,n1,n2,sigma)
x=0:0.01:10;
[l,k]=size(x);
k=k/2;
weight_init=0.1*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*sigma^2));
for i=1:sp
    thal_sp(i,i:i+2*n1)=weight_init;
end
%-------------------------------------
l4=sp-2*n2;
%----------------------------------------------------
%----------------------------------------
%connections from thalamic to layer 4 (flat)
for i=1:l4
    for j=1:2*n1+1
       thal_l4(i,j+i-1)=0.02;
    end
%     thal_l4(i,i:i+8)=2*weight_init;
end
%-------------------------------------------
 weight_init=0.5+0.5*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*((2*n1+2)/4)^2));
 for i=1:l4
     sp_l4(i,i:i+2*n2)=0.11;
 end

 
    


        
        
    
