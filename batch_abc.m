dir1='C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tun_ssa';
tms=[1000 2000 3000 4000 5000];
freq=[8 7 6 5 4 3 2];
for iii=1:length(tms)
    for jjj=1:length(freq)
        cd(dir1);
        eval(sprintf('mkdir %s_9-%s',num2str(tms(iii)),num2str(freq(jjj))));
        result_file=strcat(dir1,'\',num2str(tms(iii)),'_9-',num2str(freq(jjj)));
        Tir=tms(iii);
        f2=freq(jjj);
        cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model');
        ssa_abc;
        fprintf('%s_9-%s complete',num2str(tms(iii)),num2str(freq(jjj)));
        clearvars -except dir1 tms freq iii jjj
    end
end
flo=dir;
figure
for iii=3:37
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tun_ssa');
    subplot(5,7,iii-2)
    cd(flo(iii,1).name);
    list=dir;
    a={list.name};
    a(1:2)=[]; a(end)=[];
    str=sprintf('%s#',a{:});
    num=sscanf(str, '%*6c_%*u-%*3c-%*u_%*4c_%u.mat#');
    [dummy, index] = sort(num);
    b = a(index);
    std_spk=zeros(1,6500/50);
    load(list(end).name);
    for i=1:length(b)
        load(b{1,i});
        c=l4_inh_fire_time(f1-8,:);
        %d=l4_fire_time(f2-4,:);
        c=c(c>0);
        %d=d(d>0);
        for j=1:length(c);
            spks1(i,j)=c(j);
            std_spk(fix(c(j)/(token_len))+1)=std_spk(fix(c(j)/(token_len))+1)+1;
        end
        %     for j=1:length(d);
        %         spks2(i,j)=d(j);
        %         dev_spk(fix(d(j)/(token_len+iti))+1)=dev_spk(fix(d(j)/(token_len+iti))+1)+1;
        %     end
    end
    load(list(length(list),1).name);
    std_spk=std_spk./2;
    plot(std_spk);
    clearvars -except iii flo
end