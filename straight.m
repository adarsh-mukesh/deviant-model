function [str_mat]=straight(old_mat)
las=find(not(old_mat(1,:)));
for i=1:size(old_mat,1)
    j=1;
    for k=1:las(1,1)-1
        str_mat(i,j)=old_mat(i,i+k-1);
        j=j+1;
    end
end
end