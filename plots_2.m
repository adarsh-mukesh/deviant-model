for i=1:size(thal_fire_time,1)
    a=thal_fire_time(i,:);
    a=a(a>0);
    g_thal(i,:)=stim(10,2,0,20,0.001,a);
end

for i=1:size(sp_fire_time,1)
    a=sp_fire_time(i,:);
    a=a(a>0);
    g_sp(i,:)=stim(10,2,0,20,0.001,a);
end

figure
wt1=thal_l4_inh(21,21:29,20000);
wt1=wt1(wt1>0);
wt2=sp_l4_inh(21,21:29,20000);
wt2=wt2(wt2>0);
%raster(a);
plot(v_l4_inh(21,1:5000),'k');
hold on
plot((wt1*g_thal(25:33,1:5000)),'b');
hold on
plot((wt2*g_sp(21:29,1:5000)),'r');
hold on
plot(0.21*ones(1,5000),'g');

a=thal_fire_time(29,:);
a=a(a>0);
thspk(1,a)=1;
thspk(1,20000)=0;
b=sp_fire_time(25,:);
b=b(b>0);
spspk(1,b)=1;
spspk(1,20000)=0;
c=l4_fire_time(21,:);
c=c(c>0);
l4spk(1,c)=1;
l4spk(1,20000)=0;
figure
plot(xcorr(spspk,l4spk,50));
hold on
plot(xcorr(thspk,l4spk,50),'r');
amp_strength=0.015;
amp_weak=0.021;
tau_strength=10;
tau_weak=30;
wt=0.02;
thlfire=0;
l4lfire=0;
for i=1:max(a(end),c(end))
    if ismember(i,a)
        thlfire=i;
        if l4lfire~=0
            wt=wt*(1-amp_weak*exp(-(i-l4lfire)/tau_weak));
        end
    end
    if ismember(i,c)
        l4lfire=i;
        if thlfire~=0
            wt=wt*(1+amp_strength*exp(-(i-thlfire)/tau_strength));
        end
    end
end
wt
co='bgrky';
a=[10 17 24 31];
figure
for i=1:length(a)
    plot(thal_sp(a(i),:,1),co(i));
    hold on
end
for i=1:length(wt_thal_l4)
    a(i)=wt_thal_l4(1,i).wt(5);
end
plot(a,'k');
% tuning check
tun_dir='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\testing';
f=dir(tun_dir);
f(1:2)=[]; f(end)=[];
a={f.name};
str=sprintf('%s#',a{:});
num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
[dummy, index] = sort(num);
b = a(index);
spt=[]; l4it=[]; l4et=[];
for i=1:length(b);
    load(b{1,i});
    spst=sp_fire_time(9,:);
    l4inhst=l4_inh_fire_time(5,:);
    l4exst=l4_ex_fire_time;
    spt=[spt length(spst(spst>0))];
    l4it=[l4it length(l4inhst(l4inhst>0))];
    l4et=[l4et length(l4exst(l4exst>0))];
end
figure
plot(spt,'r')
hold on
plot(l4it,'b');
hold on
plot(l4et,'k');

figure
s=[];
for i=1:size(thal_fire_time,1)
    a=thal_fire_time(i,:);
    s=[s length(a(a>0))];
end
plot(s);
figure
s=[];
for i=1:size(sp_fire_time,1)
    a=sp_fire_time(i,:);
    s=[s length(a(a>0))];
end
plot(s);
figure
s=[];
for i=1:size(l4_inh_fire_time,1)
    a=l4_inh_fire_time(i,:);
    s=[s length(a(a>0))];
end
plot(s);




    
    

        






