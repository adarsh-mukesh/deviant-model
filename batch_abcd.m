dir1='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tun_ssa_2';
tms=[1000 2000 3000 4000 5000];
freq=[12 11 10 9 8 7 6];
for iii=1:length(tms)
    for jjj=1:length(freq)
        cd(dir1);
        eval(sprintf('mkdir %s_13-%s',num2str(tms(iii)),num2str(freq(jjj))));
        result_file=strcat(dir1,'\',num2str(tms(iii)),'_13-',num2str(freq(jjj)));
        Tir=tms(iii);
        f2=freq(jjj);
        cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model');
        ssa_abcd;
        fprintf('%s_13-%s complete\n',num2str(tms(iii)),num2str(freq(jjj)));
        clearvars -except dir1 tms freq iii jjj
    end
end

flo=dir;
st=cell(35,1);
de=cell(35,1);
fit_pts=cell(35,2);
figure
for fn=3:length(flo)
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tun_ssa_abcd_no_inh');
    subplot(5,7,fn-2);
    cd(flo(fn,1).name);
    list=dir;
    a={list.name};
    a(1:2)=[]; a(end)=[];
    str=sprintf('%s#',a{:});
    num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
    [dummy, index] = sort(num);
    b = a(index);
    std_spk=zeros(1,6500/50);
    dev_spk=zeros(1,6500/50);
    load(list(end).name);
    for fm=1:length(b)
        load(b{1,fm});
        c=l4_ex_fire_time(f1-12,:);
        d=sp_fire_time(f1-4,:);
        c=c(c>0);
        d=d(d>0);
        for jn=1:length(c);
            spks1(fm,jn)=c(jn);
            if c(jn)~=6500
                std_spk(fix(c(jn)/(token_len))+1)=std_spk(fix(c(jn)/(token_len))+1)+1;
            else
                std_spk(fix(c(jn)/(token_len)))=std_spk(fix(c(jn)/(token_len)))+1;
            end
                
        end
        for jn=1:length(d);
            spks2(fm,jn)=d(jn);
            if d(jn)~=6500
                dev_spk(fix(d(jn)/(token_len))+1)=dev_spk(fix(d(jn)/(token_len))+1)+1;
            else
                dev_spk(fix(d(jn)/(token_len)))=dev_spk(fix(d(jn)/(token_len)))+1;
            end
        end
    end
    load(list(length(list),1).name);
    std_spk=std_spk./2;
    st{fn-2,1}=std_spk;
    dev_spk=dev_spk./2;
    de{fn-2,1}=dev_spk;
    x1=[1:6:31]';
    fit_pts{fn-2,1}=(std_spk(x1))';
    x2=[1:6:31]';
    fit_pts{fn-2,2}=(dev_spk(x2))';
    plot(dev_spk,'r');
    hold on 
    plot(std_spk);
    clearvars -except fn flo st de fit_pts
end
coef=zeros(35,6,2);
figure
decay_const=cell(35,1);
for fn=1:size(fit_pts,1)
    subplot(5,7,fn);
    f1=fit([0:300:1500]',fit_pts{fn,1},'exp1');
    con1=confint(f1);
    f2=fit([0:300:1500]',fit_pts{fn,2},'exp1');
    con2=confint(f2);
    coef(fn,1,1)=f1.a; coef(fn,2:3,1)=(con1(:,1))';
    coef(fn,4,1)=f1.b; coef(fn,5:6,1)=(con1(:,2))';
    coef(fn,1,2)=f2.a; coef(fn,2:3,2)=(con2(:,1))';
    coef(fn,4,2)=f2.b; coef(fn,5:6,2)=(con2(:,2))';
    plot([1:1:6],coef(fn,1,1)*exp([0:300:1500].*coef(fn,4,1)),'b');
    t=text(5,30,num2str(-1/coef(fn,4,1)));
    set(t,'FontSize',8); set(t,'Color','b');
    hold on
    plot([1:1:6],coef(fn,1,2)*exp([0:300:1500].*coef(fn,4,2)),'r');
    k=text(5,35,num2str(-1/coef(fn,4,2)));
    set(k,'FontSize',8); set(k,'Color','r');
    decay_const{fn,1}=[-1/coef(fn,4,1) -1/coef(fn,4,1)];
end
coeff=cell(5,4);
for i=1:5
    for j=1:4
        if isempty(str_grp_devpos{i,j})~=1
        for k=1:length(str_grp_devpos{i,j})
            dp=str_grp_devpos{i,j}{1,k};
            f=fit((str_grp_tmpos{i,j}{1,k}(1:dp-1))',(str_grp_ssatc{i,j}{1,k}{1,1}(1:dp-1))','exp1');
            coeff{i,j}{1,k}(1,1:2)=[f.a f.b];
            coeff{i,j}{1,k}(2:3,1:2)=confint(f);
            dc{i,j}(1,k)=-1/f.b;
        end
        end
    end
end
a=dc{1,2}; a=a(a>0);
b=dc{1,4}; b=b(b>0);
c=dc{5,2}; c=c(c>0);
d=dc{5,4}; d=d(d>0);
figure
subplot(2,2,1)
hist(a,100);
subplot(2,2,2)
hist(b,100);
subplot(2,2,3)
hist(c,100);
subplot(2,2,4)
hist(d,100);





    
