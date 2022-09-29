sp_all=[];l4_inh_all=[];l4_ex_all=[];
cd('J:\adarsh_model\results\csi-trial_new');
ff=dir;
ff(1:2)=[];
n={ff.name};
indx=cellfun(@(v)str2num(v(12:end)),n);
[a,ind]=sort(indx);
ff=ff(ind,:);
for hh=1:9
cd(ff(hh).name)
fl=dir;fl(1:2)=[];
load(fl(end).name)
clearvars -except f1 f2 n1 n2 n3 iti token_len step fl bin_sz jump sp t_Sim step sp_all l4_inh_all l4_ex_all ff hh
fl(end)=[];
n={fl.name};
if hh==10
    indx=cellfun(@(v)str2num(v(30:end-4)),n);
else
    indx=cellfun(@(v)str2num(v(29:end-4)),n);
end
[a,ind]=sort(indx);
fl=fl(ind,:);
n_toks=t_Sim/step;
kappa=1;
sp_std_spk=zeros(2,sp); sp_dev_spk=zeros(2,sp);
l4_inh_std_spk=zeros(2,sp-2*n2); l4_inh_dev_spk=zeros(2,sp-2*n2);
l4_ex_std_spk=zeros(2,sp-2*(n3+n2)); l4_ex_dev_spk=zeros(2,sp-2*(n3+n2));

for ii=1:n_toks
    for jj=1:2
        load(fl((length(fl)/2)*(jj-1)+ii).name);
        frpt=fr(:,1:(iti+token_len):size(fr,2));
        max_fr=max(frpt(:,1));
        dev_pos=find(fr(f2*(2-jj)+f1*(jj-1),:)==max_fr);
        std_pos=find(fr(f1*(2-jj)+f2*(jj-1),:)==max_fr);
        num_std_tok=length(std_pos)/token_len;
        num_dev_tok=length(dev_pos)/token_len;
        for kk=1:sp
            frr=sp_fire_time(kk,:);
            frr=frr(frr>0);
            sp_std_spk(jj,kk)=sp_std_spk(jj,kk)+length(find(ismember(frr,std_pos)==1));
            sp_dev_spk(jj,kk)=sp_dev_spk(jj,kk)+length(find(ismember(frr,dev_pos)==1));
        end
        for kk=1:sp-2*n2
            frr=l4_inh_fire_time(kk,:);
            frr=frr(frr>0);
            l4_inh_std_spk(jj,kk)=l4_inh_std_spk(jj,kk)+length(find(ismember(frr,std_pos)==1));
            l4_inh_dev_spk(jj,kk)=l4_inh_dev_spk(jj,kk)+length(find(ismember(frr,dev_pos)==1));
        end
        for kk=1:sp-2*(n2+n3)
            frr=l4_inh_fire_time(kk,:);
            frr=frr(frr>0);
            l4_ex_std_spk(jj,kk)=l4_ex_std_spk(jj,kk)+length(find(ismember(frr,std_pos)==1));
            l4_ex_dev_spk(jj,kk)=l4_ex_dev_spk(jj,kk)+length(find(ismember(frr,dev_pos)==1));
        end
    end
end
csi_sp=(sum(sp_dev_spk./num_dev_tok)-sum(sp_std_spk./num_std_tok))./(sum(sp_dev_spk./num_dev_tok)+sum(sp_std_spk./num_std_tok));
csi_l4_inh=(sum(l4_inh_dev_spk./num_dev_tok)-sum(l4_inh_std_spk./num_std_tok))./(sum(l4_inh_dev_spk./num_dev_tok)+sum(l4_inh_std_spk./num_std_tok));
csi_l4_ex=(sum(l4_ex_dev_spk./num_dev_tok)-sum(l4_ex_std_spk./num_std_tok))./(sum(l4_ex_dev_spk./num_dev_tok)+sum(l4_ex_std_spk./num_std_tok));
sp_all=[sp_all;csi_sp];
l4_inh_all=[l4_inh_all;csi_l4_inh];
l4_ex_all=[l4_ex_all;csi_l4_ex];
cd ..
end  

figure
dff=3;
plot(sp_all(:,17-dff),'b')
hold on
plot(l4_inh_all(:,13-dff),'r')
hold on
plot(l4_ex_all(:,9-dff),'k')




    
       
            
        
        
    
    