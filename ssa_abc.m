%declare_parameters--------------------------------
global sp n1 n2
sp=9;
n1=4;
n2=4;
sigma=4;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
%stim_param----------
f1=9;
%f2=5;
token_len=50;%milliseconds
iti=250;%milliseconds
num_token=15;
dev_pos=7;
tun_wid=1;
gap=2000;%milliseconds
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
v_th_l4_inh=0.06;
u_th=1;
u_sp=1;
v_spike=10;
beta=5;
ref_length=20;
Tre=0.9;%in milliseconds
Tei=5.3;%in milliseconds
%Tir=2000;%in milliseconds
%--------------------------------------------------
%weight_generation---------------------------------
l4_inh=sp-2*n2;
thal=sp+2*n1;
thal_l4_inh=zeros(l4_inh,l4_inh+2*n1);
sp_l4_inh=zeros(l4_inh,sp);
[thal_sp, thal_l4_inh, sp_l4_inh]=weight_gen_abc(sigma);% weight matrix
%thal_l4_inh(:,:,1)=a;
%--------------------------------------------------------------
thal_fire_time=[];
sp_fire_time=zeros(sp,1);
l4_inh_fire_time=zeros(l4_inh,1);
iter_num=1;
x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,1,1)=1;
x_sp=zeros(sp,3,step*(1/dt)); x_sp(:,1,1)=1;
for step_num=1:step:t_Sim 
    %fr=max(20*exp(-(([1:57]-kappa).^2)/(2*0.5^2)),2);
    x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,1,1)=1;
    x_sp=zeros(sp,3,step*(1/dt)); x_sp(:,1,1)=1;
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
    %g_thal=ones(thal,step*(1/dt)));
    thal_fire_time(:,~any(thal_fire_time,1))=[];
        
    %-------------------------------------------------------------------
    %run thal to sp-----------------------------------------------------
    k=ones(sp,1);
    sp_last_spk=zeros(sp,1);
    v_sp=zeros(sp,step*(1/dt));
    g_sp=zeros(sp,step*(1/dt));
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
                del=1;
            else
                del=0;
            end
            x_sp(i,:,t)=depressing_synapses(x_sp(i,:,t-1),Tre,Tei,Tir,del);
        end
    end
%     for i=1:sp
%         [spont]=spk_gen_poss(fs,tin,tout,dt);
%         a=sp_fire_time(i,:);
%         a=a(a>0);
%         a=sort([a,spont]);
%         sp_fire_time(i,:)=0;
%         for j=1:length(a)
%             sp_fire_time(i,j)=a(j);
%         end
%     end
%     a=[];
    %fprintf('thal and sp firing complete\n');
    for i=1:sp
        [g_sp(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,sp_fire_time(i,:),u_sp);
    end
    %g_sp=ones(sp,step*(1/dt)));
    sp_fire_time(:,~any(sp_fire_time,1))=[];
    %-------------------------------------------------------------------
    k=ones(l4_inh,1);
    l4_inh_last_spk=zeros(l4_inh,1);% index of the l4_inh_fire_time matrix for each l4_inh neuron
    v_l4_inh=zeros(l4_inh,step*(1/dt));
    l4_inh_fire_time=zeros(l4_inh,500);
    num=0;
    gum=0;
    lum=0;
    kum=1;
    for t=2:step*(1/dt)
        for i=1:l4_inh
%---------------------------------------------------------------------------------------------------------------
            if l4_inh_last_spk(i,1)==0
                v_l4_inh(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i:i+2*n2,t-1));
            else
                if t-l4_inh_last_spk(i,1)<=ref_length
                    v_l4_inh(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i:i+2*n2,t-1))-beta*exp(-(t-l4_inh_last_spk(i,1)-1)/t_refr);
                else
                    v_l4_inh(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i:i+2*n2,t-1));
                end
            end
            if v_l4_inh(i,t)>=v_th_l4_inh
                l4_inh_last_spk(i,1)=t;
                l4_inh_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end  
        end
%         if size(l4_inh_fire_time,1)~=l4_inh
%             l4_inh_fire_time(l4_inh,1)=0;
%         end
    end
    l4_inh_fire_time(:,~any(l4_inh_fire_time,1))=[];
    for i=1:l4_inh
        [spont]=spk_gen_poss(fs,tin,tout,dt);
        a=l4_inh_fire_time(i,:);
        a=a(a>0);
        a=sort([a,spont]);
        l4_inh_fire_time(i,:)=0;
        for j=1:length(a)
            l4_inh_fire_time(i,j)=a(j);
        end
    end
    a=[];
    cd(result_file);
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(iter_num));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_inh_fire_time thal_sp sp_l4_inh thal_l4_inh v_l4_inh',resultfile))
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model');
    %fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
end
 cd(result_file);
 save sum
 cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model')
