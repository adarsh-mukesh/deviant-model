tic
%declare_parameters--------------------------------
sp=49;
n1=4;
n2=4;
sigma=4;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=4;
t_refr=2;%milliseconds (to be put in same units)
t_Sim=5*57;% in seconds ; total simulation time seconds (to be put in same units)
step=5;%time in one step (in seconds)
dt=0.001;%seconds (to be put in same units)
fr=10;%Hz (to be put in same units)  
fs=1.5;
%fs_w=0.2;
lps=0.1;% leak par sample
v_th_sp=0.11;
v_th_l4=0.2;
v_spike=10;
u_th=1;
u_sp=0.5;
%amp_strength=0.015;
amp_strength=0;
%amp_weak=0.021;
amp_weak=0;
tau_strength=13;%in millisecond
tau_weak=20;%in millisecond
beta=10;
ref_length=20;
%--------------------------------------------------
%weight_generation---------------------------------
l4=sp-2*n2;
thal=sp+2*n1;
thal_l4=zeros(l4,l4+2*n1,step/dt);
sp_l4=zeros(l4,sp,step/dt);
thal_sp=zeros(sp,thal,step/dt);
[thal_sp, thal_l4(:,:,1), sp_l4(:,:,1)]=weight_gen(sp,n1,n2,sigma);% weight matrix
%--------------------------------------------------------------
wt_thal_l4=struct;
wt_sp_l4=struct;
time_step=1;
thal_fire_time=zeros(57,1);
sp_fire_time=zeros(sp,1);
l4_fire_time=zeros(l4,1);
thal_fire_time_las=[];
sp_fire_time_las=zeros(sp,1);
l4_fire_time_las=zeros(l4,1);
iter_num=1;
stop_l4=zeros(l4,1);
kappa=1;
salk=cell(49,1);
b=ones(49,1);
for step_num=1:step:t_Sim 
    fr=max(10*exp(-(([1:57]-kappa).^2)/(2*1^2)),2);
    if step_num~=1
        thal_l4(:,:,1)=last_thal_l4;
        sp_l4(:,:,1)=last_sp_l4;
        thal_fire_time=[];
        sp_fire_time=zeros(sp,1);
        l4_fire_time=zeros(l4,1);
    end
    %stim generation-----------------------------------------------
    tin=step_num;
    tout=step_num+step;
    for i=1:thal
        [spk_mat] = spk_gen_poss(fr(i),tin,tout,dt);
        %[spk_mat] = spk_gen_poss_nonhom(fr,tin,tout,dt);
        spk_mat=sort(spk_mat);
        for j=1:size(spk_mat,2)
            thal_fire_time(i,j)=spk_mat(1,j);
        end
    end
    for i=1:thal
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);%thalamic conductances
    end
    if step_num~=1;
        thal_fire_time=[thal_fire_time_las, thal_fire_time]; % adding spikes from previous iter
    end
    %-------------------------------------------------------------------
    %run thal to sp-----------------------------------------------------
    for t=1:step/dt
        for i=1:sp
            v_sp_input(i,t)=0;
            for j=i:2*n1+i
                v_sp_input(i,t)=v_sp_input(i,t)+g_thal(j,t)*thal_sp(i,j);
            end
        end
    end
    % for t=1:tSim/dt-2
    %     for i=1:sp
    %         v_sp_decay(i,t+1)=v_sp_input(i,t+1)-0.1*v_sp_input(i,t);
    %     end
    % end
    k=ones(sp,1);
    sp_last_spk=zeros(sp,1);
    v_sp=zeros(sp,step/dt-1);
    for t=2:step/dt
        for i=1:sp
            if sp_last_spk(i,1)==0
                v_sp(i,t)=(1-lps)*v_sp_input(i,t-1);
            else
                if t-sp_last_spk(i,1)<=ref_length
                    v_sp(i,t)=(1-lps)*v_sp_input(i,t-1)-beta*exp(-(t-sp_last_spk(i,1)-1)/t_refr);
                else
                    v_sp(i,t)=(1-lps)*v_sp_input(i,t-1);
                end
            end
            if v_sp(i,t)>=v_th_sp
                sp_last_spk(i,1)=t;
                sp_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end
        end
    end
    if size(sp_fire_time,1)~=sp
        sp_fire_time(sp,1)=0;
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
    if step_num~=1
        sp_fire_time=[sp_fire_time_las, sp_fire_time]; %adding spikes from previous iter
    end
    %-------------------------------------------------------------------
    k=ones(l4,1);
    l4_last_spk=zeros(l4,1);% index of the l4_fire_time matrix for each l4 neuron
    v_l4=zeros(l4,step/dt-1);
    if step_num~=1
        k(1:l4,1)=2;
        l4_fire_time(:,1)=l4_fire_time_las; % adding spikes from previous iter
    end
    num=0;
    gum=0;
    lum=0;
    kum=1;
    for t=2:step/dt
        for i=1:l4
            
            %if l4_las_fire(i,1)~=0
            %kappa=99;
            % now check which sp and thal neurons have fired on this time
            % and see the difference from the last time when the
            % corresponding l4 fired
            if stop_l4(i,1)==0 %----------------testing for pre synaptic spike--------------------------------------------------
                for j=i+n2:i+n2+2*n1   %thal to l4 checking
                    if ismember(t,thal_fire_time(j,:))==1
                        if l4_last_spk(i,1)~=0
                            thal_l4(i,j-n2,t)=update_weak(l4_last_spk(i,1),t,thal_l4(i,j-n2,t-1),amp_weak,tau_weak);
                            gum=gum+1;
                            salk{i,1}(b(i,1),1)=thal_l4(i,j-n2,t-1);
                            salk{i,1}(b(i,1),2)=thal_l4(i,j-n2,t);
                            salk{i,1}(b(i,1),3)=t;
                            salk{i,1}(b(i,1),4)=-1;
                            b(i,1)=b(i,1)+1;
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
                        sp_l4(i,j,t)=update_weak(l4_last_spk(i,1),t,sp_l4(i,j,t-1),amp_weak,tau_weak);
                    if sp_l4(i,j,t)~=sp_l4(i,j,t-1)
                        lum=lum+1;
                    end
                    else
                        sp_l4(i,j,t)=sp_l4(i,j,t-1);
                    end
                else
                    sp_l4(i,j,t)=sp_l4(i,j,t-1);
                end
            end%---------------------------------------------------------------------------------------------------------------
            wt=thal_l4(i,i:i+2*n1,t);
            input=g_thal(i+n2:i+n2+2*n1,t);
            v_l4_input(i,t)=wt*input;
            wt=sp_l4(i,i:i+2*n2,t);
            input=g_sp(i:i+2*n2,t);
            v_l4_input(i,t)=v_l4_input(i,t)+wt*input;
            if l4_last_spk(i,1)==0
                v_l4(i,t)=(1-lps)*v_l4_input(i,t-1);
            else
                if t-l4_last_spk(i,1)<=ref_length
                    v_l4(i,t)=(1-lps)*v_l4_input(i,t-1)-beta*exp(-(t-l4_last_spk(i,1)-1)/t_refr);
                else
                    v_l4(i,t)=(1-lps)*v_l4_input(i,t-1);
                end
            end
            if v_l4(i,t)>=v_th_l4
                if stop_l4(i,1)==0
                    for j=i+n2:i+n2+2*n1
                        up_lim=0.3;
                        thal_l4(i,j-n2,t)=update_strength(thal_fire_time(j,:),t,thal_l4(i,j-n2,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                        num=num+1;
                        salk{i,1}(b(i,1),1)=thal_l4(i,j-n2,t-1);
                        salk{i,1}(b(i,1),2)=thal_l4(i,j-n2,t);
                        salk{i,1}(b(i,1),3)=t;
                        salk{i,1}(b(i,1),4)=1;
                        b(i,1)=b(i,1)+1;
                        if thal_l4(i,j-n2,t)==up_lim
                            stop_l4(i,1)=0;
                        end
                    end
                else
                    for j=i+n2:i+n2+2*n1
                        thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
                    end
                end
                for j=i:i+2*n2
                    up_lim=0.11;
                    sp_l4(i,j,t)=update_strength(sp_fire_time(j,:),t,sp_l4(i,j,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                end
                l4_last_spk(i,1)=t;
                l4_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            else
%                 for j=i+n2:i+n2+2*n1
%                     thal_l4(i,j-n2,t)=thal_l4(i,j-n2,t-1);
%                 end
%                 for j=i:i+2*n2
%                     sp_l4(i,j,t)=sp_l4(i,j,t-1);
%                 end
            end  
        end
%         if size(l4_fire_time,1)~=l4
%             l4_fire_time(l4,1)=0;
%         end
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
    last_thal_l4=thal_l4(:,:,step/dt);
    wt_thal_l4(1,time_step).wt=thal_l4(:,:,step/dt);
    last_sp_l4=sp_l4(:,:,step/dt);
    wt_sp_l4(1,time_step).wt=sp_l4(:,:,step/dt);
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
    
    for l=1:size(l4_fire_time,1)
        m=l4_fire_time(l,:);
        m=m(m~=0);
        if isempty(m)==1
            l4_fire_time_las(l,1)=0;
        else
            l4_fire_time_las(l,1)=m(end)-step/dt;
        end
        
    end
    
    cd('C:\adarsh_model\results\tuning_wala');
    str=date;
    resultfile=strcat('result_',str,'_iter_00',num2str(time_step));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_fire_time thal_sp sp_l4 thal_l4',resultfile))
    cd('C:\adarsh_model');
    time_step=time_step+1;
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
   kappa=kappa+1; 
end
toc
