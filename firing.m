function [ poisson_fire,weight ] = firing(num_sub, fam_num )
k=num_sub+(fam_num-1);
for ii=1:k
    poisson_fire(k,:)=poissonSpikeGen(30,1,1);
end
 for m=1:(fam_num+1)/2
   w(m)=1/(((fam_num+1)/2)+1-m);
   w(fam_num-m+1)=1/(((fam_num+1)/2)+1-m);
 end

w=w/sum(w);
l=0;
for i=1:num_sub
    for j=1:fam_num
        weight(i,j+l)=w(j);
    end 
    l=l+1;
end
        
        
    
%for j=1:num_sub
  %  input=poisson_fire(:,j:j+fam_num-1);
    
 %   fire_sub(j)=sum(poisson_fire(:,j:j+fam_num-1),2);
%end
        
        

end

