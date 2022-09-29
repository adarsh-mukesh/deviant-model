function[]=revercheck(mat)
mat(:,1)=[];
for i=1:size(mat,1)
    m=mat(i,:);
    m=m(m~=0);
    if issorted(m)==0
        fprintf('load hai %i isme\n',i);
    end
end
    