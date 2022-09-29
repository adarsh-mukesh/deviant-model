sp_all=[];l4_all=[];
cd('J:\adarsh_model\results\csi_bin');
ff=dir;
ff(1:2)=[];
n={ff.name};
indx=cellfun(@(v)str2num(v(12:end)),n);
[a,ind]=sort(indx);
ff=ff(ind,:);
for hh=1:4
cd(ff(hh).name)
fl=dir;fl(1:2)=[];
load(fl(end).name)
clearvars -except spont iti token_len step fl t_sim step sp_all l4_all ff hh
fl(end)=[];
n={fl.name};
if hh==10
    indx=cellfun(@(v)str2num(v(30:end-4)),n);
else
    indx=cellfun(@(v)str2num(v(27:end-4)),n);
end
[a,ind]=sort(indx);
fl=fl(ind,:);
n_toks=t_sim/step;
kappa=1;
sp_std_spk=zeros(2,1); sp_dev_spk=zeros(2,1);
l4_std_spk=zeros(2,1); l4_dev_spk=zeros(2,1);

for ii=1:n_toks
    for jj=1:2
        load(fl((length(fl)/2)*(jj-1)+ii).name);
        frpt=fr(:,1:(iti+token_len):size(fr,2));
        max_fr=max(frpt(:,1));
        dev_pos=find(fr(2*(2-jj)+1*(jj-1),:)==max_fr);
        std_pos=find(fr(1*(2-jj)+2*(jj-1),:)==max_fr);
        num_std_tok=length(std_pos)/token_len;
        num_dev_tok=length(dev_pos)/token_len;
        
            frr=sp_fire_time(1,:);
            frr=frr(frr>0);
            sp_std_spk(jj,1)=sp_std_spk(jj,1)+length(find(ismember(frr,std_pos)==1));
            sp_dev_spk(jj,1)=sp_dev_spk(jj,1)+length(find(ismember(frr,dev_pos)==1));
        
        
            frr=l4_fire_time(1,:);
            frr=frr(frr>0);
            l4_std_spk(jj,1)=l4_std_spk(jj,1)+length(find(ismember(frr,std_pos)==1));
            l4_dev_spk(jj,1)=l4_dev_spk(jj,1)+length(find(ismember(frr,dev_pos)==1));
        
    end
end
csi_sp=(sum(sp_dev_spk./num_dev_tok)-sum(sp_std_spk./num_std_tok))./(sum(sp_dev_spk./num_dev_tok)+sum(sp_std_spk./num_std_tok));
csi_l4=(sum(l4_dev_spk./num_dev_tok)-sum(l4_std_spk./num_std_tok))./(sum(l4_dev_spk./num_dev_tok)+sum(l4_std_spk./num_std_tok));
sp_all=[sp_all;csi_sp];
l4_all=[l4_all;csi_l4];
cd ..
end  