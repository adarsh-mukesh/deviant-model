dir1=pwd;
dir2='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tuning_wala_results';
cd(dir2);
f=dir;
list=dir(fullfile(dir2,'*mat'));
file_name1={list.name};
str=sprintf('%s#',file_name1{:});
num=sscanf(str,'%*6c_%*11c_%*4c_%u.mat#');
[dum idx]=sort(num);
for jj=1:length(dum)
    file_name2(jj)=file_name1(idx(jj));
end


color='mcrgbk';
for i=1:length(file_name2)
    load(file_name2{1,i});
    b=[];
    for j=1:size(l4_fire_time,1)
        a=l4_fire_time(j,:);
        b(j)=length(a(a>0));
    end
    f_t_l4(i,:)=b;
end
for i=1:length(file_name2)
    load(file_name2{1,i});
    b=[];
    for j=1:size(sp_fire_time,1)
        a=sp_fire_time(j,:);
        b(j)=length(a(a>0));
    end
    f_t_sp(i,:)=b;
end
 cd(dir1);
 figure
 for i=1:size(f_t_l4,2)-1
     c=f_t_l4(i+4:i+12,i);
     subplot(5,8,i)
     plot(c);
 end
 figure
 for i=1:size(f_t_sp,2)
     c=f_t_sp(i:i+8,i);
     subplot(7,7,i);
     plot(c);
 end
% 

