tic
%declare_parameters--------------------------------
global sp n1 n2 n3
sp=49;
n1=4;
n2=4;
n3=4;
sigma=4;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
t_Sim=20;% in seconds ; total simulation time seconds (to be put in same units)
step=20;%time in one step (in seconds)
dt=0.001;%seconds (to be put in same units)
fr=10;%Hz (to be put in same units)
fs=1.5;
%fs_w=0.2;
lps=0.1;% leak par sample
v_th_sp=0.4;
v_th_l4_inh=0.4;
v_th_l4_ex=0.4;
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
%--------------------------------------------------
%weight_generation---------------------------------
l4_inh=sp-2*n2;
l4_ex=l4_inh-2*n3;
thal=sp+2*n1;
thal_l4_inh=zeros(l4_inh,l4_inh+2*n1,step/dt);
thal_l4_ex=zeros(l4_ex,l4_ex+2*n1,step/dt);
sp_l4_inh=zeros(l4_inh,sp,step/dt);
sp_l4_ex=zeros(l4_ex,l4_ex+2*n2,step/dt);
%l4_inh_ex=zeros(l4_ex,l4_inh,step/dt);
[thal_sp, thal_l4_inh(:,:,1), sp_l4_inh(:,:,1), thal_l4_ex(:,:,1), sp_l4_ex(:,:,1), l4_inh_ex]=weight_gen_abcd(sigma);% weight matrix
%thal_l4(:,:,1)=a;
%--------------------------------------------------------------
wt_thal_l4_inh=struct;
wt_thal_l4_ex=struct;
wt_sp_l4_inh=struct;
wt_sp_l4_ex=struct;
%wt_l4_inh_ex=struct;
time_step=1;
thal_fire_time=[];
sp_fire_time=zeros(sp,1);
l4_inh_fire_time=zeros(l4_inh,1);
l4_ex_fire_time=zeros(l4_ex,1);
thal_fire_time_las=[];
sp_fire_time_las=zeros(sp,1);
l4_inh_fire_time_las=zeros(l4_inh,1);
l4_ex_fire_time_las=zeros(l4_ex,1);
iter_num=1;
stop_l4_inh=zeros(l4_inh,1);
stop_l4_ex=zeros(l4_ex,1);
kappa=1;
for step_num=1:step:t_Sim 
    %fr=max(20*exp(-(([1:57]-kappa).^2)/(2*0.5^2)),2);
    if step_num~=1
        thal_l4_inh(:,:,1)=last_thal_l4_inh;
        sp_l4_inh(:,:,1)=last_sp_l4_inh;
        thal_l4_ex(:,:,1)=last_thal_l4_ex;
        sp_l4_ex(:,:,1)=last_sp_l4_ex;
        %l4_inh_ex(:,:,1)=last_l4_inh_ex;
    end
    %stim generation-----------------------------------------------
    tin=step_num;
    tout=step_num+step;
    g_thal=zeros(thal,step/dt);
    thal_fire_time=zeros(thal,500);
    for i=1:thal
        [spk_mat] = spk_gen_poss(fr,tin,tout,dt);
        thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);%thalamic conductances
    end
    %g_thal=ones(thal,step/dt);
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
        for i=1:sp
            if sp_last_spk(i,1)==0
                v_sp(i,t)=(1-lps)*(thal_sp(i,i:i+2*n1)*g_thal(i:i+2*n1,t-1));
            else
                if t-sp_last_spk(i,1)<=ref_length
                    v_sp(i,t)=(1-lps)*(thal_sp(i,i:i+2*n1)*g_thal(i:i+2*n1,t-1))-beta*exp(-(t-sp_last_spk(i,1)-1)/t_refr);
                else
                    v_sp(i,t)=(1-lps)*(thal_sp(i,i:i+2*n1)*g_thal(i:i+2*n1,t-1));
                end
            end
            if v_sp(i,t)>=v_th_sp
                sp_last_spk(i,1)=t;
                sp_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end
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
    fprintf('thal and sp firing complete\n');
    for i=1:sp
        [g_sp(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,sp_fire_time(i,:),u_sp);
    end
    %g_sp=ones(sp,step/dt);
    if step_num~=1
        sp_fire_time=[sp_fire_time_las, sp_fire_time]; %adding spikes from previous iter
    end
    sp_fire_time(:,~any(sp_fire_time,1))=[];
    %-------------------------------------------------------------------
    k=ones(l4_inh,1);
    g_l4_inh=zeros(l4_inh,step/dt);
    l4_inh_last_spk=zeros(l4_inh,1);% index of the l4_fire_time m  atrix for each l4 neuron
    v_l4_inh=zeros(l4_inh,step/dt);
    l4_inh_fire_time=zeros(l4_inh,500);
    if step_num~=1
        k(1:l4_inh,1)=2;
        l4_inh_fire_time(:,1)=l4_inh_fire_time_las; % adding spikes from previous iter
    end
    for t=2:step/dt
        for i=1:l4_inh
            % now check which sp and thal neurons have fired on this time
            % and see the difference from the last time when the
            % corresponding l4 fired
            if stop_l4_inh(i,1)==0 %----------------testing for pre synaptic spike--------------------------------------------------
                for j=i+n2:i+n2+2*n1   %thal to l4 checking
                    if ismember(t,thal_fire_time(j,:))==1
                        if l4_inh_last_spk(i,1)~=0
                            thal_l4_inh(i,j-n2,t)=update_weak(l4_inh_last_spk(i,1),t,thal_l4_inh(i,j-n2,t-1),amp_weak,tau_weak);
                        else
                            thal_l4_inh(i,j-n2,t)=thal_l4_inh(i,j-n2,t-1);
                        end
                    else
                        thal_l4_inh(i,j-n2,t)=thal_l4_inh(i,j-n2,t-1);
                    end
                end
            else
                for j=i+n2:i+n2+2*n1
                    thal_l4_inh(i,j-n2,t)=thal_l4_inh(i,j-n2,t-1);
                end
            end
            for j=i:i+2*n2     %sp to l4 checking
                if ismember(t,sp_fire_time(j,:))==1
                    if l4_inh_last_spk(i,1)~=0
                        sp_l4_inh(i,j,t)=update_weak(l4_inh_last_spk(i,1),t,sp_l4_inh(i,j,t-1),amp_weak,tau_weak);
                    else
                        sp_l4_inh(i,j,t)=sp_l4_inh(i,j,t-1);
                    end
                else
                    sp_l4_inh(i,j,t)=sp_l4_inh(i,j,t-1);
                end
            end
%---------------------------------------------------------------------------------------------------------------
            if l4_inh_last_spk(i,1)==0
                v_l4_inh(i,t)=(1-lps)*(thal_l4_inh(i,i:i+2*n1,t-1)*g_thal(i+n2:i+n2+2*n1,t-1)+sp_l4_inh(i,i:i+2*n2,t-1)*g_sp(i:i+2*n2,t-1)); 
            else
                if t-l4_inh_last_spk(i,1)<=ref_length
                    v_l4_inh(i,t)=(1-lps)*(thal_l4_inh(i,i:i+2*n1,t-1)*g_thal(i+n2:i+n2+2*n1,t-1)+sp_l4_inh(i,i:i+2*n2,t-1)*g_sp(i:i+2*n2,t-1))-beta*exp(-(t-l4_inh_last_spk(i,1)-1)/t_refr);
                else
                    v_l4_inh(i,t)=(1-lps)*(thal_l4_inh(i,i:i+2*n1,t-1)*g_thal(i+n2:i+n2+2*n1,t-1)+sp_l4_inh(i,i:i+2*n2,t-1)*g_sp(i:i+2*n2,t-1));
                end
            end
            if v_l4_inh(i,t)>=v_th_l4_inh
                if stop_l4_inh(i,1)==0
                    for j=i+n2:i+n2+2*n1
                        up_lim=0.4;
                        thal_l4_inh(i,j-n2,t)=update_strength(thal_fire_time(j,:),t,thal_l4_inh(i,j-n2,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                        if thal_l4_inh(i,j-n2,t)==up_lim
                            stop_l4_inh(i,1)=0;
                        end
                    end
                else
                    for j=i+n2:i+n2+2*n1
                        thal_l4_inh(i,j-n2,t)=thal_l4_inh(i,j-n2,t-1);
                    end
                end
                for j=i:i+2*n2
                    up_lim=0.11;
                    sp_l4_inh(i,j,t)=update_strength(sp_fire_time(j,:),t,sp_l4_inh(i,j,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                end
                l4_inh_last_spk(i,1)=t;
                l4_inh_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end  
        end
    end
    l4_inh_fire_time(:,~any(l4_inh_fire_time,1))=[];
    for i=1:l4_inh
        if step_num==1
            [g_l4_inh(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,l4_inh_fire_time(i,:),u_l4);
        else
            [g_l4_inh(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,l4_inh_fire_time(i,2:end),u_l4);
        end
    end
    k=ones(l4_ex,1);
    l4_ex_last_spk=zeros(l4_ex,1);% index of the l4_fire_time m  atrix for each l4 neuron
    v_l4_ex=zeros(l4_ex,step/dt);
    l4_ex_fire_time=zeros(l4_ex,500);
    if step_num~=1
        k(1:l4_ex,1)=2;
        l4_ex_fire_time(:,1)=l4_ex_fire_time_las; % adding spikes from previous iter
    end
    for t=2:step/dt
        for i=1:l4_ex
            for j=i+n2+n3:i+n2+n3+2*n1   %thal to l4 checking
                if ismember(t,thal_fire_time(j,:))==1
                    if l4_ex_last_spk(i,1)~=0
                        thal_l4_ex(i,j-n2-n3,t)=update_weak(l4_ex_last_spk(i,1),t,thal_l4_ex(i,j-n2-n3,t-1),amp_weak,tau_weak);
                    else
                        thal_l4_ex(i,j-n2-n3,t)=thal_l4_ex(i,j-n2-n3,t-1);
                    end
                else
                    thal_l4_ex(i,j-n2-n3,t)=thal_l4_ex(i,j-n2-n3,t-1);
                end
            end
            for j=i+n3:i+n3+2*n2     %sp to l4 checking
                if ismember(t,sp_fire_time(j,:))==1
                    if l4_ex_last_spk(i,1)~=0
                        sp_l4_ex(i,j-n3,t)=update_weak(l4_ex_last_spk(i,1),t,sp_l4_ex(i,j-n3,t-1),amp_weak,tau_weak);
                    else
                        sp_l4_ex(i,j-n3,t)=sp_l4_ex(i,j-n3,t-1);
                    end
                else
                    sp_l4_ex(i,j-n3,t)=sp_l4_ex(i,j-n3,t-1);
                end
            end
%             for j=i:i+2*n3
%                 if ismember(t,l4_inh_fire_time(j,:))==1
%                     if l4_ex_last_spk(i,:)~=0
%                         l4_inh_ex(i,j,t)=update_weak(l4_ex_last_spk(i,1),t,l4_inh_ex(i,j,t-1),amp_weak,tau_weak);
%                     else
%                         l4_inh_ex(i,j,t)=l4_inh_ex(i,j,t-1);
%                     end
%                 else
%                     l4_inh_ex(i,j,t)=l4_inh_ex(i,j,t-1);
%                 end
%             end
            if l4_ex_last_spk(i,1)==0
                v_l4_ex(i,t)=(1-lps)*(thal_l4_ex(i,i:i+2*n1,t-1)*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+sp_l4_ex(i,i:i+2*n2,t-1)*g_sp(i+n3:i+n3+2*n2,t-1)-l4_inh_ex(i,i:i+2*n3)*g_l4_inh(i:i+2*n3,t-1));
            else
                if t-l4_ex_last_spk(i,1)<=ref_length
                    v_l4_ex(i,t)=(1-lps)*(thal_l4_ex(i,i:i+2*n1,t-1)*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+sp_l4_ex(i,i:i+2*n2,t-1)*g_sp(i+n3:i+n3+2*n2,t-1)-l4_inh_ex(i,i:i+2*n3)*g_l4_inh(i:i+2*n3,t-1))-beta*exp(-(t-l4_ex_last_spk(i,1)-1)/t_refr);
                else
                    v_l4_ex(i,t)=(1-lps)*(thal_l4_ex(i,i:i+2*n1,t-1)*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+sp_l4_ex(i,i:i+2*n2,t-1)*g_sp(i+n3:i+n3+2*n2,t-1)-l4_inh_ex(i,i:i+2*n3)*g_l4_inh(i:i+2*n3,t-1));
                end
            end
            if v_l4_ex(i,t)>=v_th_l4_ex
                if stop_l4_ex(i,1)==0
                    for j=i+n2+n3:i+n2+n3+2*n1
                        up_lim=0.4;
                        thal_l4_ex(i,j-n2-n3,t)=update_strength(thal_fire_time(j,:),t,thal_l4_ex(i,j-n2-n3,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                        if thal_l4_ex(i,j-n2-n3,t)==up_lim
                            stop_l4_ex(i,1)=0;
                        end
                    end
                else
                    for j=i+n2+n3:i+n2+n3+2*n1
                        thal_l4_ex(i,j-n2-n3,t)=thal_l4_ex(i,j-n2-n3,t-1);
                    end
                end
                for j=i+n3:i+n3+2*n2
                    up_lim=0.11;
                    sp_l4_ex(i,j-n3,t)=update_strength(sp_fire_time(j,:),t,sp_l4_ex(i,j-n3,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                end
%                 for j=i:i+2*n3
%                     up_lim=0.1;
%                     l4_inh_ex(i,j,t)=update_strength(l4_inh_fire_time(j,:),t,l4_inh_ex(i,j,t-1),amp_strength,amp_weak,tau_strength,up_lim);
%                 end
                l4_ex_last_spk(i,1)=t;
                l4_ex_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end
        end
    end
        
%     for i=1:l4
%         [spont]=spk_gen_poss(fs,tion,tout,dt);
%         a=l4_fire_time(i,:);
%         a=a(a>0);
%         a=sort([a,spont]);
%         l4_fire_time(i,:)=0;
%         for j=1:length(a)
%             l4_fire_time(i,j)=a(j);
%         end
%     end
%     a=[];
    last_thal_l4_inh=thal_l4_inh(:,:,step/dt);
    wt_thal_l4_inh(1,time_step).wt=thal_l4_inh(:,:,step/dt);
    last_sp_l4_inh=sp_l4_inh(:,:,step/dt);
    wt_sp_l4_inh(1,time_step).wt=sp_l4_inh(:,:,step/dt);
    last_thal_l4_ex=thal_l4_ex(:,:,step/dt);
    wt_thal_l4_ex(1,time_step).wt=thal_l4_ex(:,:,step/dt);
    last_sp_l4_ex=sp_l4_ex(:,:,step/dt);
    wt_sp_l4_ex(1,time_step).wt=sp_l4_ex(:,:,step/dt);
    %last_l4_inh_ex=l4_inh_ex(:,:,step/dt);
%     wt_l4_inh_ex(1,time_step).wt=l4_inh_ex(:,:,step/dt);
    thal_fire_time_las=zeros(thal,1);
    for l=1:size(thal_fire_time,1)
        m=thal_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)~=1
            thal_fire_time_las(l,1)=m(end)-step/dt;
            
        end
    end
    for l=1:size(sp_fire_time,1)
        m=sp_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)~=1
            sp_fire_time_las(l,1)=m(end)-step/dt;
        end
    end
    
    for l=1:size(l4_inh_fire_time,1)
        m=l4_inh_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)==1
            l4_inh_fire_time_las(l,1)=0;
        else
            l4_inh_fire_time_las(l,1)=m(end)-step/dt;
        end
    end
    
    for l=1:size(l4_ex_fire_time,1)
        m=l4_ex_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)==1
            l4_ex_fire_time_las(l,1)=0;
        else
            l4_ex_fire_time_las(l,1)=m(end)-step/dt;
        end
    end
    cd('C:\adarsh_model\results\tuning_wala');
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(time_step));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_inh_fire_time l4_ex_fire_time thal_sp sp_l4_inh sp_l4_ex thal_l4_inh thal_l4_ex l4_inh_ex v_l4_inh v_l4_ex',resultfile))
    cd('C:\adarsh_model');
    
    time_step=time_step+1;
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
   kappa=kappa+1;
end
cd('C:\adarsh_model\results\tuning_wala');
save sum
cd('C:\adarsh_model');
toc
