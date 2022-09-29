function [ fire_time,weight ] = ini_fire(num_sub, fam_num )
num=num_sub+(fam_num-1);
for ii=1:num
    poisson_fire=poissonSpikeGen(20,1,1,0.001);
    y=find(poisson_fire);
    [m,n]=size(y);
    for jj=1:n
        fire_time(ii,jj)=y(jj);
    end
    
    
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
%for ii=1:num_sub
 %   fire_1=zeros(1,1000);
  %  fire_fam=fire_time(ii:fam_num,:);
   % [m,n]=size(fire_fam);
    %index=ones(1,m);
   % for uu=1:1000
   % while fire_fam(:,index(1):index(m))~=0
     %       for cz=1:m
      %          if fire_fam(cz,index(cz))==uu
       %             fire_1(fire_fam(cz,index(cz)))=fire_1(fire_fam(cz,index(cz)))+1;
        %            index(cz)=index(cz)+1;
                   

           
           
    
    
    
        
        

end

