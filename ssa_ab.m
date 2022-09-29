tic
%declare_parameters--------------------------------
global sp n1
sp=1;
n1=4;
sigma=4;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
%stim_param----------
f1=9;
f2=5;
token_len=50;%milliseconds
iti=250;%milliseconds
num_token=15;
dev_pos=7;
tun_wid=1;
gap=5000;%milliseconds
%------------------
step=(15*(token_len+iti)+gap)/1000;%time in one step (in seconds)
t_Sim=40*step;% in seconds ; total simulation time seconds (to be put in same units)
dt=0.001;%seconds (to be put in same units)
fr=fire_ssa_gen(f1,f2,token_len,iti,num_token,dev_pos,tun_wid,gap);
%fr=pure_ssatone(f1,f2,token_len,iti,num_token,dev_pos);
fs=1;
%fs_w=0.2;
lps=0.1;% leak par sample
v_th_sp=0.04;
u_th=1;
u_sp=1;
v_spike=10;
beta=5;
ref_length=20;
Tre=0.9;%in milliseconds
Tei=5.3;%in milliseconds
Tir=5000;%in milliseconds
%--------------------------------------------------
%weight_generation---------------------------------
thal=sp+2*n1;
thal_sp=zeros(sp,sp+2*n1);
[thal_sp]=weight_gen_ab(sigma);% weight matrix
%thal_l4(:,:,1)=a;
%------------------------------------------------------------
thal_fire_time=[];
sp_fire_time=zeros(sp,1);
sp_fire_time_las=zeros(sp,1);
iter_num=1;
for step_num=1:step:t_Sim 
    %fr=max(20*exp(-(([1:57]-kappa).^2)/(2*0.5^2)),2);
    x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,1,1)=1;
    %stim generation-----------------------------------------------
    tin=step_num;
    tout=step_num+step;
    g_thal=zeros(thal,step*(1/dt));
    thal_fire_time=zeros(thal,500);
    for i=1:thal
        [spk_mat] = spk_gen_poss_nonhom(fr(i,:),tin,tout,dt);
        thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);%thalamic conductances
    end
    thal_fire_time(:,~any(thal_fire_time,1))=[];
    %-------------------------------------------------------------------
    %run thal to sp-----------------------------------------------------
    k=ones(sp,1);
    sp_last_spk=zeros(sp,1);
    v_sp=zeros(sp,step*(1/dt));
    sp_fire_time=zeros(sp,500);
    for t=2:step*(1/dt)
        for i=1:thal
            if ismember(t,thal_fire_time(i,:))
                del=1;
            else
                del=0;
            end
            x_thal(i,:,t)=depressing_synapses(x_thal(i,:,t-1),Tre,Tei,Tir,del);
        end
        for i=1:sp
            if sp_last_spk(i,1)==0
                v_sp(i,t)=(1-lps)*(thal_sp(i,i:i+2*n1)*g_thal(i:i+2*n1,t-1));
            else
                if t-sp_last_spk(i,1)<=ref_length
                    v_sp(i,t)=(1-lps)*((squeeze(x_thal(i:i+2*n1,2,t-1))'.*thal_sp(i,i:i+2*n1))*g_thal(i:i+2*n1,t-1))-beta*exp(-(t-sp_last_spk(i,1)-1)/t_refr);
                else
                    v_sp(i,t)=(1-lps)*((squeeze(x_thal(i:i+2*n1,2,t-1))'.*thal_sp(i,i:i+2*n1))*g_thal(i:i+2*n1,t-1));
                end
            end
            if v_sp(i,t)>=v_th_sp
                sp_last_spk(i,1)=t;
                sp_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end
        end
    end
    for i=1:sp
        [spont]=spk_gen_poss(fs,tin,tout,dt);
        a=sp_fire_time(i,:);
        a=a(a>0);
        a=sort([a,spont]);
        sp_fire_time(i,:)=0;
        for j=1:length(a)
            sp_fire_time(i,j)=a(j);
        end
    end
    a=[];
    fprintf('thal and sp firing complete\n');
    %-------------------------------------------------------------------
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tuning_wala_results\testing');
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(iter_num));
    eval(sprintf('save %s thal_fire_time sp_fire_time v_sp',resultfile))
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model');
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
end
 cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running\tuning_wala_results\testing');
 save sum
toc
