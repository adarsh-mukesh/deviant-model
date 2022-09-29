% batch-----------------------------------
figure
for bn=9:14
    subplot(1,6,bn-8)
    k=[];
    for vn=1:length(wt_sp_l4_inh)
        k=[k ;wt_sp_l4_inh(1,vn).wt(bn,bn:bn+8)];
    end
    imagesc(k)
end
%-----------------------------------------
% iter------------------------------------
figure
for bn=1:40
    subplot(5,8,bn)
    k=squeeze(sp_l4_inh(bn,bn+4,:));
    plot(k);
end
%-----------------------------------------
%
figure
for bn=1:40
    subplot(5,8,bn)
    k=wt_thal_l4(1,length(wt_thal_l4)).wt(bn,bn:bn+8);
    plot(k)
end
figure
for bn=17:32
    subplot(4,4,bn-16)
    k=wt_thal_l4(1,120).wt(bn,bn:bn+16)
    plot(k)
end
% l4 fire rate------------------------------
figure
b=[];
for i=1:41
    a=l4_inh_fire_time(i,:);
    b(i)=length(a(a>0));
end
plot([zeros(1,8) b zeros(1,8)]);
%--------------------------------------------
hold on
b=[];
for i=1:49
    a=sp_fire_time(i,:);
    b(i)=length(a(a>0));
end
plot([zeros(1,4) b zeros(1,4)],'k');
%-------------------------------------------
hold on
a=[];
for i=1:33
    b=l4_ex_fire_time(i,:);
    a(i)=length(b(b>0));
end
plot([zeros(1,12) a zeros(1,12)],'r');
% batch--------------------------------------
figure
for i=1:40
    subplot(5,8,i);
    plot(squeeze(sp_l4(i,i+4,:)));
end

figure
n1=4;
n2=4;
sigma=3;
for i=21:21
    %subplot(5,8,i);
    a=thal_l4_inh(i,i:i+8,20000);
    plot(a);
%     hold on
%     plot(0.2*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*sigma^2)),'r');
%     hold on
%     plot(mean(a)*ones(1,9),'b');
%     hold on
%     plot(mean(0.2*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*sigma^2)))*ones(1,9),'r');
end


figure
for bn=1:16
    subplot(5,8,bn)
    k=[];
    for vn=1:length(wt_thal_l4)
        k=[k wt_thal_l4(1,vn).wt(bn,bn+4)];
    end
    plot(k);
end
% ts_weak
ts=[];
m=1;
for i=1:l4
    a=l4_fire_time(i,:);
    a=a(a~=0);
    for j=i:i+2*n2
        b=sp_fire_time(j,:);
        b=b(b~=0);
        for k=1:length(b)
            l=a(a<b(k));
            if isempty(l)
                ts(m,k)=-0.5;
            else
                ts(m,k)=b(k)-l(end);
            end
        end
        m=m+1;
    end
end

testing=cell(49,2);
for i=1:49
    a=sp_fire_time(i,:);
    a=a(a>0);
    k=0;
    for j=1:length(a)
        b=v_sp_input(a(j)-1);
        b=b-[0.0 0.1 0.0]*g_thal(i+3:i+5,a(j)-1);
        if b<v_th_sp
            k=k+1;
        end
        testing{i,1}=k/length(a);
    end
end
testing=cell(41,2);
for i=1:41
    a=l4_fire_time(i,:);
    a=a(a>0);
    k=[];
    for j=1:length(a)
        for ii=-1:1:1
            b=sp_fire_time(i+n2+ii,:);
            b=b(b>0);
            c=b(b<a(j));
            if isempty(c)~=1
                testing{i,1}(ii+2,j)=a(j)-c(end);
            else
                testing{i,1}(ii+2,j)=-1;
            end
        end
    end
end
for i=1:41
    a=l4_fire_time(i,:);
    a=a(a>0);
    k=[];
    for j=1:length(a)
        for ii=-1:1:1
            b=sp_fire_time(i+n2+ii,:);
            b=b(b>0);
            c=b(b>a(j));
            if isempty(c)~=1
                testing{i,2}(ii+2,j)=c(1)-a(j);
            else
                testing{i,2}(ii+2,j)=-1;
            end
        end
    end
end
figure
for ii=1:16
    k=[];
    subplot(4,4,ii)
    for i=1:length(wt_thal_l4_ex)
        a=wt_thal_l4_ex(1,i).wt(ii,:);
        a=a(a~=0);
        k(i,:)=a;
    end
    imagesc(k);
end
figure
for i=1:9;
    subplot(3,3,i)
    k=[];
    for j=1:180
        k=[k wt_thal_l4(1,j).wt(21,20+i)];
    end
    plot(k);
end
k=[];
for i=1:360
    k=[k wt_thal_l4(1,i).wt(1,5)];
end

figure
plot(k);
list=dir('K:\adarsh_model\results');
a={list.name};
a(1:2)=[]; a(end)=[];
str=sprintf('%s#',a{:});
num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
[dummy, index] = sort(num);
b = a(index);
std_spk=zeros(1,18000/50);
dev_spk=zeros(1,18000/50);
cd('K:\adarsh_model\results');
load(list(end).name);
spks1=zeros(length(a),40);
spks2=zeros(length(a),40);
ss=[];
tt=[];
for i=1:length(b)
    load(b{1,i});
    c=l4_inh_fire_time(f1-8,:);
    d=sp_fire_time(f1-4,:);
    c=c(c>0);
    d=d(d>0);
    for j=1:length(c);
        if c(j)~=6500
            spks1(i,j)=c(j);
            std_spk(fix(c(j)/(token_len))+1)=std_spk(fix(c(j)/(token_len))+1)+1;
        end
    end
    for j=1:length(d);
        if d(j)~=6500
            spks2(i,j)=d(j);
            dev_spk(fix(d(j)/(token_len))+1)=dev_spk(fix(d(j)/(token_len))+1)+1;
        end
    end
    s=[];
    for i=1:size(thal_fire_time,1)
        a=thal_fire_time(i,:);
        s=[s length(a(a>0))];
    end
    ss=[ss;s];
     t=[];
    for i=1:size(sp_fire_time,1)
        a=sp_fire_time(i,:);
        t=[t length(a(a>0))];
    end
    tt=[tt;t];
    %     v_lf=zeros(1,6500);
    %     %x_thal1=ones(25,3,6500);
    %     for t=2:6500
    %         v_lf(t)=(1-lps)*((squeeze(x_thal(1+n2+n3:1+2*n1+n2+n3,2,t-1))'.*thal_l4_inh(1,1:1+2*n1))*g_thal(1+n2+n3:1+n2+n3+2*n1,t-1)+(squeeze(x_sp(1+n3:1+n3+2*n2,2,t-1))'.*sp_l4_inh(1,1:1+2*n2))*g_sp(1+n3:1+n3+2*n2,t-1)-(squeeze(x_l4_inh(1:1+2*n3,2,t-1))'.*l4_inh_ex(1,1:1+2*n3))*g_l4_inh(1:1+2*n3,t-1));
    %     end
    %     figure
    %     plot(v_lf);
end

load(list(length(list),1).name);
cd('K:\adarsh_model');
%dev_spk=length(b);
std_spk=std_spk./5;
dev_spk=dev_spk./5;
figure
 plot(std_spk);
 hold on
plot(dev_spk,'r');
figure
for i=1:num_token
    if i==dev_pos
        ed_color='r';
        f_color='r';
    else
        ed_color='c';
        f_color='c';
    end
    rectangle('position',[(iti+token_len)*(i-1) 0 token_len length(b)],'EdgeColor',ed_color,'FaceColor',f_color);
    hold on
end
raster(spks1);

%%%%%%%%%%%%%%%mean and std time constant
ac=[];bc=[];
for jj=1:4
    for ii=1:5
        int_1=str_grp_tcf{iis,jj};
        ac{ii,jj}(1,1)=nanmean(abs(int_1(:,1))); ac{ii,jj}(1,2)=nanmean(abs(int_1(:,2)));
        bc{ii,jj}(1,1)=nanstd(abs(int_1(:,1)),1,1)/sqrt(length(int_1(:,1))); bc{ii,jj}(1,2)=nanstd(abs(int_1(:,2)),1,1)/sqrt(length(int_1(:,2)));
    end
end

for pp=1:2 %%% 1=TN 'r' 2=NT 'b'
    for ii=1:4
        str_tc=[];
        for jj=1:5
            dd=[ac{jj,ii}(1,pp) ;bc{jj,ii}(1,pp)];
            str_tc=[str_tc dd] ;
        end
        subplot(4,1,ii); hold on
        if pp==1
            errorbar(str_tc(1,:),str_tc(2,:),'r')
        elseif pp==2
            errorbar(str_tc(1,:),str_tc(2,:),'k')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
k=1;
data_path='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\madhur_system\re_17_tuning';
list = dir('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\madhur_system\re_17_tuning');
a={list.name};
a(1:2)=[];a(end)=[];
str=sprintf('%s#',a{:});
num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
[dummy, index] = sort(num);
b = a(index);
d=[];
cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\madhur_system\re_17_tuning');
for i=5:length(b)-4
    load(b{1,i});
    subplot(7,7,k);
    for j=1:size(sp_fire_time,1)
        a=sp_fire_time(j,:);
        a=a(a>0);
        c(j)=length(a);
    end
    a=sp_fire_time(i-4,:);
    a=a(a>0);
    d=[d length(a)];
    plot(c);
    hold on
    k=k+1;
end
cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model');

data_path='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\madhur_system\re_18';
list = dir('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\madhur_system\re_18');
a={list.name};
a(1:2)=[];
a(end)=[];
str=sprintf('%s#',a{:});
num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
[dummy, index] = sort(num);
b = a(index);
d=[];
for i=1:length(b)
    load(b{1,i});
    d=[d; thal_l4(:,:,20000)];
    ft=l4_fire_time(21,:);
    e(i)=length(ft(ft>0));
end
color='mcrgbk'
x=1:180;
figure
plotyy(x,e,x,k);
hold on
for j=1:6
    plot([j*10 j*10],[0 300],color(j));
end

figure
for j=1:6
    plot(wt_thal_l4(1,j*10).wt(21,21:29),color(j));
    hold on
end

list=dir('J:\adarsh_model\results\all_freq');
a={list.name};
a(1:2)=[]; 
str=sprintf('%s#',a{:});
num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
[dummy, index] = sort(num);
b = a(index);

frt=[];
frr=[];
frtd=[];
for i=1:length(b)
    load(b{1,i});
    %frr=[frr;fr];
        c=l4_ex_fire_time(9,:);
        c=c(c>0);
        frt=[frt length(c)];
%         d=thal_fire_time(i,:);
%         d=d(d>0);
%         frtd=[frtd length(d)];
        
    
end
frt=frt./10;
i=9;
v_in=zeros(1,6500);
x_thal1=ones(25,3,6500);
for t=2:6500
    v_in(t)=(1-lps)*((squeeze(x_thal(i:i+2*n1,2,t-1))'.*thal_sp(i,i:i+2*n1))*g_thal(i:i+2*n1,t-1));
end
i=1;
v_lf=zeros(1,4500);
for t=2:4500
    v_lf(t)=(1-lps)*((squeeze(x_thal(i+n2+n3:i+2*n1+n2+n3,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+(squeeze(x_sp(i+n3:i+n3+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i+n3:i+n3+2*n2,t-1)-(squeeze(x_l4_inh(i:i+2*n3,2,t-1))'.*l4_inh_ex(i,i:i+2*n3))*g_l4_inh(i:i+2*n3,t-1));
end 
figure
for nn=1:9
    subplot(3,3,nn)
    plot(squeeze(x_thal(8+nn,2,:)));
end
pst_l4=zeros(1,18000/50);
a=l4_ex_fire_time;
a=a(a>0);
for i=1:length(a)
    pst_l4(fix(a(i)/(token_len))+1)=pst_l4(fix(a(i)/(token_len))+1)+1;
end
pst_l4=pst_l4.*20;
figure
ak=fr(f2,:);
ak=ak(1:50:end);
ak=find(ak==10);
for i=1:length(ak)
    plot([ak(i) ak(i)],[min(pst_l4) max(pst_l4)],'g');
    hold on
end
plot(pst_l4)
al=zeros(25,500);
for i=1:25
    [test]=spk_gen_poss_nonhom(fr(i,:),0,6.5,0.001);
    al(i,1:length(test))=test;
end
al(:,~any(al,1))=[];
    

    


    





    
    
    













