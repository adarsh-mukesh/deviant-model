tic
%declare_parameters--------------------------------
global sp n1 n2
sp=17;
n1=8;
n2=8;
sigma=2;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
% stim params-------------------------------
f1=17;
f2=12;
token_len=50;%milliseconds
iti=250;%milliseconds
num_token=15;
dev_pos=7;
tun_wid=2;
%gap=2000;%milliseconds
%-------------------------------------------
step=60*(iti+token_len)/1000;
%step=(15*(token_len+iti)+gap)/1000;%time in one step (in seconds)
t_Sim=2*step;% in seconds ; total simulation time seconds (to be put in same units)
dt=0.001;%seconds (to be put in same units)
%fr=fire_ssa_gen(f1,f2,token_len,iti,num_token,dev_pos,tun_wid,gap);%Hz (to be put in same units)
%fr=fire_oddball_seq(f1,f2,token_len,iti,tun_wid,step,dt);
fs=0.5;
%fs_w=0.2;
lps=0.1;% leak par sample
v_th_sp=0.12;
v_th_l4=0.12;
u_th=1;
u_sp=1;
u_l4=1;
v_spike=10;
amp_strength=0.015;
%amp_strength=0;
amp_weak=0.021;
%amp_weak=0;
tau_strength=13;%in millisecond
tau_weak=20;%in millisecond
beta=5;
ref_length=20;
Tre_thal=0.9; Tre_sp=0.9;%in milliseconds
Tei_thal=10; Tei_sp=27;%in milliseconds
Tir_thal=3000; Tir_sp=3000;%in milliseconds
%--------------------------------------------------
%weight_generation---------------------------------
l4=sp-2*n2;
thal=sp+2*n1;
thal_sp=zeros(sp,thal);
thal_l4=zeros(l4,l4+2*n1,step*(1/dt));
sp_l4=zeros(l4,sp,step*(1/dt));
[thal_sp, thal_l4(:,:,1), sp_l4(:,:,1)]=weight_gen_abc(sigma);% weight matrix
%thal_l4(:,:,1)=a;
%--------------------------------------------------------------
thal_fire_time=[];
sp_fire_time=zeros(sp,1);
l4_fire_time=zeros(l4,1);
x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,1,1)=1;
x_sp=zeros(sp,3,step*(1/dt)); x_sp(:,1,1)=1;
wt_thal_l4=struct;
wt_sp_l4=struct;
stop_l4=zeros(l4,1);
iter_num=1;
tir=Tir_thal;
rs=zeros(l4,1);
for step_num=1%:step:t_Sim
    fr=fire_oddball_seq(f1,f2,token_len,iti,tun_wid,step,dt);
    %fr=50*exp(-(([1:25]-iter_num).^2)/(2*tun_wid^2));
   %fr=fire_allfreq_seq(token_len,iti,tun_wid,step,dt);
    if step_num~=1
        x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,:,1)=last_x_thal;
        x_sp=zeros(sp,3,step*(1/dt)); x_sp(:,:,1)=last_x_sp;
        sp_l4=zeros(l4,sp,step*(1/dt)); sp_l4(:,:,1)=last_sp_l4;
        thal_l4=zeros(l4,l4+2*n1,step*(1/dt)); thal_l4(:,:,1)=last_thal_l4;
    end
    %stim generation-----------------------------------------------
    tin=step_num;
    tout=step_num+step;
    g_thal=zeros(thal,step/dt);
    thal_fire_time=zeros(thal,500);
    for i=1:thal
        [spk_mat] = spk_gen_poss_nonhom(fr(i,:),tin,tout,dt);
        %[spk_mat]=spk_gen_poss(fr(i),tin,tout,dt);
        thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);%thalamic conductances
    end
    if step_num~=1;
        thal_fire_time=[thal_fire_time_las, thal_fire_time]; % adding spikes from previous iter
    end
    thal_fire_time(:,~any(thal_fire_time,1))=[];
    
    %-------------------------------------------------------------------
    %run thal to sp-----------------------------------------------------
    k=ones(sp,1);
    sp_last_spk=zeros(sp,1);
    v_sp=zeros(sp,step/dt);
    g_sp=zeros(sp,step/dt);
    sp_fire_time=zeros(sp,500);
    for t=2:step/dt
        for i=1:thal
            if ismember(t,thal_fire_time(i,:))
                del=1;
            else
                del=0;
            end
            x_thal(i,:,t)=depressing_synapses(x_thal(i,:,t-1),Tre_thal,Tei_thal,Tir_thal,del);
        end
        for i=1:sp
            if sp_last_spk(i,1)==0
                v_sp(i,t)=(1-lps)*((squeeze(x_thal(i:i+2*n1,2,t-1))'.*thal_sp(i,i:i+2*n1))*g_thal(i:i+2*n1,t-1));
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
            x_sp(i,:,t)=depressing_synapses(x_sp(i,:,t-1),Tre_sp,Tei_sp,Tir_sp,del);
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
    %fprintf('thal and sp firing complete\n');
    for i=1:sp
        [g_sp(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,sp_fire_time(i,:),u_sp);
    end
    if step_num~=1
        sp_fire_time=[sp_fire_time_las, sp_fire_time]; %adding spikes from previous iter
    end
    sp_fire_time(:,~any(sp_fire_time,1))=[];
    %-------------------------------------------------------------------
    k=ones(l4,1);
    g_l4=zeros(l4,step/dt);
    l4_last_spk=zeros(l4,1);% index of the l4_fire_time m  atrix for each l4 neuron
    v_l4=zeros(l4,step/dt);
    l4_fire_time=zeros(l4,500);
    if iter_num~=1
        k=2*ones(l4,1);
        l4_fire_time(:,1)=l4_fire_time_las; 
    end
    for t=2:step/dt
        for i=1:l4
            % now check which sp and thal neurons have fired on this time
            % and see the difference from the last time when the
            % corresponding l4 fired
            if stop_l4(i,1)==0 %----------------testing for pre synaptic spike--------------------------------------------------
                for j=i+n2:i+n2+2*n1   %thal to l4 checking
                    if ismember(t,thal_fire_time(j,:))==1
                        if l4_last_spk(i,1)~=0
                            thal_l4(i,j-n2,t)=update_weak(l4_last_spk(i,1),t,thal_l4(i,j-n2,t-1),amp_weak,tau_weak);
                        else
                            thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
                        end
                    else
                        thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
                    end
                end
            else
                for j=i+n2:i+n2+2*n1
                    thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
                end
            end
            for j=i:i+2*n2     %sp to l4 checking
                if ismember(t,sp_fire_time(j,:))==1
                    if l4_last_spk(i,1)~=0
                        %sp_l4(i,j,t)=update_weak(l4_last_spk(i,1),t,sp_l4(i,j,t-1),amp_weak,tau_weak);
                        sp_l4(i,j,t)=sp_l4(i,j,t-1);
                    else
                        sp_l4(i,j,t)=sp_l4(i,j,t-1);
                    end
                else
                    sp_l4(i,j,t)=sp_l4(i,j,t-1);
                end
            end
            %---------------------------------------------------------------------------------------------------------------
            if l4_last_spk(i,1)==0
                v_l4(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4(i,i:i+2*n1,t-1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4(i,i:i+2*n2,t-1))*g_sp(i:i+2*n2,t-1));
            else
                if t-l4_last_spk(i,1)<=ref_length
                    v_l4(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4(i,i:i+2*n1,t-1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4(i,i:i+2*n2,t-1))*g_sp(i:i+2*n2,t-1))-beta*exp(-(t-l4_last_spk(i,1)-1)/t_refr);
                else
                    v_l4(i,t)=(1-lps)*((squeeze(x_thal(i+n2:i+2*n1+n2,2,t-1))'.*thal_l4(i,i:i+2*n1,t-1))*g_thal(i+n2:i+n2+2*n1,t-1)+(squeeze(x_sp(i:i+2*n2,2,t-1))'.*sp_l4(i,i:i+2*n2,t-1))*g_sp(i:i+2*n2,t-1));
                end
            end
            if v_l4(i,t)>=v_th_l4
                if stop_l4(i,1)==0
                    for j=i+n2:i+n2+2*n1
                        up_lim=0.6;
                        thal_l4(i,j-n2,t)=update_strength(thal_fire_time(j,:),t,thal_l4(i,j-n2,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                        if thal_l4(i,j-n2,t)==up_lim
                            stop_l4(i,1)=0;
                        end
                    end
                    a=thal_l4(i,i:i+2*n1,t);
                    a=[a sp_l4(i,i:i+2*n2,t)];
                    if ismember(up_lim,a)
                        rs(i,1)=1;
                    end
                    if rs(i,1)==1
                        b=thal_l4(i,i:i+2*n1,t);
                        thal_l4(i,i:i+2*n1,t)=b./1; %norm(a);
                        b=sp_l4(i,i:i+2*n2,t);
                        sp_l4(i,i:i+2*n2,t)=b./1;%norm(a);
                    end
                else
                    for j=i+n2:i+n2+2*n1
                        thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
                    end
                end
                for j=i:i+2*n2
                    up_lim=0.11;
                    %sp_l4(i,j,t)=update_strength(sp_fire_time(j,:),t,sp_l4(i,j,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                    sp_l4(i,j,t)=sp_l4(i,j,t-1);
                end
                l4_last_spk(i,1)=t;
                l4_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            
                
            end
        end
    end
    for i=1:l4
        [spont]=spk_gen_poss(fs,tin,tout,dt);
        a=l4_fire_time(i,:);
        a=a(a>0);
        a=sort([a,spont]);
        l4_fire_time(i,:)=0;
        for j=1:length(a)
            l4_fire_time(i,j)=a(j);
        end
    end
    a=[];
    l4_fire_time(:,~any(l4_fire_time,1))=[];
    %lf_fr_rate=[lf_fr_rate length(a)];
    %a=[];
    last_thal_l4=thal_l4(:,:,end);
    wt_thal_l4(1,iter_num).wt=thal_l4(:,:,end);
    last_sp_l4=sp_l4(:,:,end);
    wt_sp_l4(1,iter_num).wt=sp_l4(:,:,end);
    last_x_thal=x_thal(:,:,end);
    last_x_sp=x_sp(:,:,end);
    thal_fire_time_las=zeros(thal,1);
    for l=1:size(thal_fire_time,1)
        m=thal_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)~=1
            thal_fire_time_las(l,1)=m(end)-step/dt;
        end
    end
    sp_fire_time_las=sp_last_spk-step/dt;
    l4_fire_time_las=l4_last_spk-step/dt;
    cd('J:\adarsh_model\results');
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(iter_num));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_fire_time thal_sp sp_l4 thal_l4 g_thal g_sp g_l4 v_l4 v_sp x_thal x_sp fr',resultfile))
    cd('J:\adarsh_model');
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
    tir=tir*0.999
end
cd('J:\adarsh_model\results');
save sum
cd('J:\adarsh_model');
toc