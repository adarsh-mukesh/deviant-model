tic
clear
d1=pwd;
time_stamp=500;
cd('J:\adarsh_model\results\bin_std_no_rise');
load('J:\adarsh_model\results\bin_std_no_rise\sum.mat')
f1=dir;
f1(1:2)=[];f1(end)=[];
n={f1.name};
indx=cellfun(@(v)str2num(v(25:end-4)),n);
[a,ind]=sort(indx);
f1=f1(ind,:);
load(f1(time_stamp).name);
clearvars -except d1 time_stamp tir thal_sp wt_thal_l4 wt_sp_l4
cd(d1);
result_file='J:\adarsh_model\results\csi_trial_no_std\time_stamp_100';
iti=250;
token_len=50;
dt=0.001;
fs=0.5;
spont_win=500;%milliseconds
num_tok=15;
dev_loc=7;
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
%step=60*(iti+token_len)/1000;% seconds
step=(spont_win+num_tok*(token_len+iti))/1000;
u_th=1;
u_sp=1;
beta=5;
ref_length=20;
t_sim=100*step;
v_th_sp=0.05;
v_th_l4=0.08;
x_thal=zeros(2,3,step*(1/dt)); x_thal(:,1,1)=1;
x_sp=zeros(step*(1/dt),3); x_sp(1,1)=1;
thal_fire_time=zeros(2,500);
sp_fire_time=zeros(1,500);
l4_fire_time=zeros(1,500);
stop_l4=0;
%init weights
thal_sp=0.2*ones(2,1);
thal_l4=0.22*ones(2,1);%wt_thal_l4(:,time_stamp);
sp_l4=0.04;%0.11;%wt_sp_l4(:,time_stamp);
tir_thal=100;
tir_sp=6000;
iter_num=1;
lps=0.1;
ts=0;
Tei_thal=10;
Tei_sp=27;
plast=zeros(2,1);
for pair=1:2
    for step_num=1:step:t_sim
        %fr=binary_fire(step,token_len,iti,dt);
        fr=binary_slice_fire(spont_win,token_len,iti,num_tok,dev_loc);
        if pair==2
            fr=fr([2 1],:);
        end
        thal_fire_time=zeros(2,500);
        x_thal=zeros(2,3,step*(1/dt)); x_thal(:,1,1)=1;
        x_sp=zeros(step*(1/dt),3); x_sp(1,1)=1;
        tin=step_num;
        tout=step_num+step;
        for i=1:2
            [spk_mat]=spk_gen_poss_nonhom(fr(i,:),tin,tout,dt);
            thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
            [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);
        end
        thal_fire_time(:,~any(thal_fire_time,1))=[];
        k_sp=1;
        k_l4=1;
        sp_last_spk=0;
        v_sp=zeros(1,step/dt);
        g_sp=zeros(1,step/dt);
        g_l4=zeros(1,step/dt);
        l4_last_spk=0;% index of the l4_fire_time m  atrix for each l4 neuron
        v_l4=zeros(1,step/dt);
        sp_fire_time=zeros(1,500);
        l4_fire_time=zeros(1,500);
        for t=2:step*(1/dt)
            for i=1:2
                if ismember(t,thal_fire_time(i,:))
                    del=1;
                else
                    del=0;
                end
                x_thal(i,:,t)=depressing_synapses(x_thal(i,:,t-1),0.9,Tei_thal,tir_thal,del);
            end
            if sp_last_spk==0
                v_sp(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_sp)'*g_thal(:,t-1));
            else
                if t-sp_last_spk<=ref_length
                    v_sp(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_sp)'*g_thal(:,t-1))-beta*exp(-(t-sp_last_spk-1)/t_refr);
                else
                    v_sp(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_sp)'*g_thal(:,t-1));
                end
            end
            if ((v_sp(1,t)>=v_th_sp) || (rand<=1-exp(-1*fs*dt)))
                sp_last_spk=t;
                sp_fire_time(1,k_sp)=t;
                k_sp=k_sp+1;
                del=1;
            else
                del=0;
            end
            x_sp(t,:)=depressing_synapses(x_sp(t-1,:),0.9,Tei_sp,tir_thal,del);
            a=sp_fire_time;
            a=a(a>t-(10*tau_syn));
            a=a(a<t);
            a=a(a>0);
            if isempty(a)
                g_sp(1,t)=0;
            else
                g_sp(1,t)=0;
                for nn=1:length(a)
                    g_sp(1,t)=g_sp(1,t)+exp(-(t-a(nn))/tau_syn);
                end
            end
            %l4 wala star
            if l4_last_spk==0
                v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,1)*g_sp(1,1));
            else
                if t-l4_last_spk<=ref_length
                    v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,1)*g_sp(1,t-1))-beta*exp(-(t-l4_last_spk-1)/t_refr);
                else
                    v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,1)*g_sp(1,t-1));
                end
            end
            if ((v_l4(1,t)>v_th_l4) || (rand<=1-exp(-1*fs*dt)))
                l4_last_spk=t;
                l4_fire_time(1,k_l4)=t;
                k_l4=k_l4+1;
                del=1;
            end
        end
        cd(result_file);
        str=date;
        resultfile=strcat('result_',num2str(time_stamp),'_',str,'_iter_',num2str(iter_num));
        eval(sprintf('save %s thal_fire_time sp_fire_time l4_fire_time sp_l4 thal_l4 g_thal g_sp v_l4 v_sp x_thal x_sp fr ',resultfile));
        cd(d1);
        fprintf('iter %i complete\n',iter_num);
        iter_num=iter_num+1;
    end
end
cd(result_file);
save sum
cd(d1)
toc