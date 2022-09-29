bin_sz=100;
jump=10;
cd('J:\adarsh_model\results\17l4_sd_rescale');
fl=dir;fl(1:2)=[];
load(fl(end).name)
clearvars -except f1 f2 n1 n2 n3 iti token_len step fl bin_sz jump sp
sp_ind=17;
l4_inh_ind=13;
l4_exp_ind=9;
fl(end)=[];
n={fl.name};
indx=cellfun(@(v)str2num(v(25:end-4)),n);
[a,ind]=sort(indx);
fl=fl(ind,:);
kappa=1;
for ii=1:jump:length(fl)-bin_sz+1
    sp_std_spk=zeros(1,sp); sp_dev_spk=zeros(1,sp);
    sp_num_std=zeros(1,sp); sp_num_dev=zeros(1,sp);
    l4_inh_std_spk=zeros(1,sp-2*n2); l4_inh_dev_spk=zeros(1,sp-2*n2);
    l4_inh_num_std=zeros(1,sp-2*n2); l4_inh_num_dev=zeros(1,sp-2*n2);
    l4_ex_std_spk=zeros(1,sp-2*(n3+n2)); l4_ex_dev_spk=zeros(1,sp-2*(n3+n2));
    l4_ex_num_std=zeros(1,sp-2*(n3+n2)); l4_ex_num_dev=zeros(1,sp-2*(n3+n2));
    for jj=ii:ii+bin_sz-1
        load(fl(jj).name);
        frpt=fr(:,1:(iti+token_len):size(fr,2));
        max_fr=max(frpt(:,1));
        dev_pos=find(fr(f2,:)==max_fr);
        std_pos=find(fr(f1,:)==max_fr);
        for kk=1:sp
            frr=sp_fire_time(kk,:);
            frr=frr(frr>0);
            sp_std_spk(1,kk)=sp_std_spk(1,kk)+length(find(ismember(frr,std_pos)==1));
            sp_dev_spk(1,kk)=sp_dev_spk(1,kk)+length(find(ismember(frr,dev_pos)==1));
            sp_num_std(1,kk)=sp_num_std(1,kk)+(length(std_pos)/token_len);
            sp_num_dev(1,kk)=sp_num_dev(1,kk)+(length(dev_pos)/token_len);
        end
        for kk=1:sp-2*n2
            frr=l4_inh_fire_time(kk,:);
            frr=frr(frr>0);
            l4_inh_std_spk(1,kk)=l4_inh_std_spk(1,kk)+length(find(ismember(frr,std_pos)==1));
            l4_inh_dev_spk(1,kk)=l4_inh_dev_spk(1,kk)+length(find(ismember(frr,dev_pos)==1));
            l4_inh_num_std(1,kk)=l4_inh_num_std(1,kk)+(length(std_pos)/token_len);
            l4_inh_num_dev(1,kk)=l4_inh_num_dev(1,kk)+(length(dev_pos)/token_len);
        end
        for kk=1:sp-2*(n2+n3)
            frr=l4_inh_fire_time(kk,:);
            frr=frr(frr>0);
            l4_ex_std_spk(1,kk)=l4_ex_std_spk(1,kk)+length(find(ismember(frr,std_pos)==1));
            l4_ex_dev_spk(1,kk)=l4_ex_dev_spk(1,kk)+length(find(ismember(frr,dev_pos)==1));
            l4_ex_num_std(1,kk)=l4_ex_num_std(1,kk)+(length(std_pos)/token_len);
            l4_ex_num_dev(1,kk)=l4_ex_num_dev(1,kk)+(length(dev_pos)/token_len);
        end
    end
    ssa_ind_sp(kappa,:)=((sp_dev_spk./sp_num_dev)-(sp_std_spk./sp_num_std))./((sp_dev_spk./sp_num_dev)+(sp_std_spk./sp_num_std));
    ssa_ind_l4_inh(kappa,:)=((l4_inh_dev_spk./l4_inh_num_dev)-(l4_inh_std_spk./l4_inh_num_std))./((l4_inh_dev_spk./l4_inh_num_dev)+(l4_inh_std_spk./l4_inh_num_std));
    ssa_ind_l4_ex(kappa,:)=((l4_ex_dev_spk./l4_ex_num_dev)-(l4_ex_std_spk./l4_ex_num_std))./((l4_ex_dev_spk./l4_ex_num_dev)+(l4_ex_std_spk./l4_ex_num_std));    
    kappa=kappa+1;
end

            
            
            
        
        
    
    