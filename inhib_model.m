%declare_parameters--------------------------------
sp=32;
n1=8;
n2=8;
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_Sim=3600;% in seconds ; total simulation time seconds (to be put in same units)
step=20;%time in one step (in seconds)
dt=0.001;%seconds (to be put in same units)
fr=10;%Hz (to be put in same units)
fs=2;
fs_w=2.2;
v_th_sp=1;
v_th_l4=2;
v_spike=10;
amp_strength=0.021;
amp_weak=0.003;
tau_strength=5;%in millisecond
tau_weak=20;%in millisecond
%--------------------------------------------------
%weight_generation---------------------------------
l4=sp-2*n2;
l4_ext=l4;
l4_inh=l4;
thal=sp+2*n1;
thal_l4=zeros(l4,l4+2*n1,step/dt);
sp_l4=zeros(l4,sp,step/dt);
thal_sp=zeros(sp,thal,step/dt);
[thal_sp, thal_l4(:,:,1), sp_l4(:,:,1)]=weight_gen(sp,n1,n2);% weight matrix
%--------------------------------------------------------------
wt_thal_l4=struct;
wt_sp_l4=struct;
time_step=1;
thal_fire_time_las=[];
sp_fire_time_las=[];
l4_fire_time_las=[];
iter_num=1;

for step_num=1:step:t_Sim
    if step_num~=1  %time_steps start here
        thal_l4_inh(:,:,1)=thal_l4_inh(:,:,step/dt);
        sp_l4_inh(:,:,1)=sp_l4_inh(:,:,step/dt);
        thal_l4_ext(:,:,1)=thal_l4_ext(:,:,step/dt);
        sp_l4_ext(:,:,1)=sp_l4_ext(:,:,step/dt);
        inh_ext(:,:,1)=inh_ext(:,:,step/dt);
        thal_fire_time=[];
        thal_fire_time(:,1)=last_thal;
        sp_fire_time=[];
        sp_fire_time(:,1)=last_sp;
        l4_inh_fire_time=[];
        l4_inh_fire_time(:,1)=last_l4_inh;
        l4_ext_fire_time=[];
        l4_ext_fire_time(:,1)=last_l4_ext;
    end
    for i=1:thal  %% thalamic inputs generated here
        [spk_mat]=spk_gen_poss(fr,tin,tout,dt);
        spk_mat=sort(spk_mat);
        for j=1:length(spk_mat)
            thal_fire_time(i,j)=spk_mat(j);
        end
    end
    if step_num~=1
        thal_fire_time=[thal_fire_time_las, thal_fire_time];
    end
    thal_last_spk=zeros(thal,1);
    for t=1:step/dt-1
        for i=1:thal
            if thal_last_spk(i)==0
                g_thal(i,t)=0;
            elseif t-thal_last_spk(i)>=t_delay
                g_thal(i,t)=exp(-(t-thal_last_spk)/tau_syn);
            end
        end
        k_sp=2*ones(sp,1);
        for i=1:sp
            v_sp_input(i,t)=0;
            for j=i:2*n1+i
                v_sp_input(i,t)=v_sp_input(i,t)+g_thal(j,t)*thal_sp(i,j);
            end
            if t==1
                if v_sp_input(i,t)>=v_sp_th
                    sp_fire_time(i,k_sp(i,1))=t;
                    v_sp(i,t)=v_spike;
                    k_sp(i,1)=k_sp(i,1)+1;
                else
                    v_sp(i,t)=v_sp_input(i,t);
                end
            else
                if t-1==sp_fire_time(i,k_sp(i,1))
                    v_sp(i,t)=0;
                else
                    if v_sp_input(i,t)>=v_sp_th
                        sp_fire_time(i,k_sp(i,1))=t;
                        v_sp(i,t)=v_spike-0.1*v_sp_input(i,t-1);
                        k_sp(i,1)=k_sp(i,1)+1;
                    else
                        v_sp(i,t)=v_sp_input(i,t)-0.1*v_sp_input(i,t-1);
                    end
                end
            end      
        end
        k_l4_inh=2*ones(l4,1);
        for i=1:l4_inh
            for j=i+n2:i+n2+2*n1
                if ismember(t,thal_fire_time(j,:))
                    last=l4_inh_fire_time(i,:);
                    last=last(last~=0);
                    thal_l4_inh(i,j-n2,t+1)=update_weak(last(end),t,thal_l4_inh(i,j-n2,t),amp_weak,tau_weak);
                else
                    thal_l4_inh(i,j-n2,t+1)=thal_l4_inh(i,j-n2,t);
                end
            end
            for j=i:i+2*n2
                if ismember(t,sp_fire_time(j,:))
                    last=l4_inh_fire_time(j,:);
                    last=last(last~=0);
                    sp_l4_inh(i,j,t+1)=update_weak(last(end),t,sp_l4_inh(i,j,t),amp_weak,tau_weak);
                else
                    sp_l4_inh(i,j,t+1)=sp_l4_inh(i,j,t);
                end
            end
            v_l4_inh_input(i,t)=0;
            for j=i:2*n2+i
                v_l4_inh_input(i,t)=v_l4_inh_input(i,t)+sp_l4_inh(i,j)*v_sp(j,t); 
            end
            for j=i:2*n1+i
                v_l4_inh_input(i,t)=v_l4_inh_input(i,t)+ thal_l4_inh(i,j)*g_thal(j+n2,t);
            end
            if t==1
                if v_l4_inh_input(i,t)>=v_l4_th
                   l4_inh_fire_time(i,k_l4_inh(i,1))=t;
                   v_l4_inh(i,t)=v_spike;
                   for j=i+n2:i+n2+2*n1
                       last=thal_fire_time(j,:);
                       last=last(last~=0);
                       last=last(last<=t);
                       thal_l4_inh(i,j,t+1)=update_strength(last(end),t,thal_l4_inh(i,j,t),amp_strength,tau_strength);
                   end
                   for j=i:i+2*n2
                       last=sp_fire_time(j,:);
                       last=last(last~=0);
                       last=last(last<=t);
                       sp_l4_inh(i,j,t+1)=update_strength(last(end),t,sp_l4_inh(i,j,t),amp_strength,tau_strength);
                   end
                   k_l4_inh(i,1)=k_l4_inh(i,1)+1;
                else
                    v_l4_inh(i,t)=v_l4_inh_input(i,t);
                end
            else
                if t-1==l4_inh_fire_time(i,k_l4_inh(i,1))
                    v_l4_inh(i,t)=0;
                else
                    if v_l4_inh_input(i,t)>=v_l4_th
                        l4_inh_fire_time(i,k_l4_inh(i,1))=t;
                        v_l4_inh(i,t)=v_spike;
                        for j=i+n2:i+n2+2*n1
                            last=thal_fire_time(j,:);
                            last=last(last~=0);
                            last=last(last<=t);
                            thal_l4_inh(i,j,t+1)=update_strength(last(end),t,thal_l4_inh(i,j,t),amp_strength,tau_strength);
                        end
                        for j=i:i+2*n2
                            last=sp_fire_time(j,:);
                            last=last(last~=0);
                            last=last(last<=t);
                            sp_l4_inh(i,j,t+1)=update_strength(last(end),t,sp_l4_inh(i,j,t),amp_strength,tau_strength);
                        end
                        k_l4_inh(i,1)=k_l4_inh(i,1)+1;
                    else
                        v_l4_inh(i,t)=v_l4_inh_input(i,t)-0.1*v_l4_inh_input(i,t-1);
                    end
                end
            end
        end
        k_l4_ext=2*ones(l4,1);
        for i=1:l4_ext
            for j=i+n2:i+n2+2*n1
                if ismember(t,thal_fire_time(j,:))
                    last=l4_ext_fire_time(i,:);
                    last=last(last~=0);
                    thal_l4_ext(i,j-n2,t+1)=update_weak(last(end),t,thal_l4_ext(i,j-n2,t),amp_weak,tau_weak);
                else
                    thal_l4_ext(i,j-n2,t+1)=thal_l4_ext(i,j-n2,t);
                end
            end
            for j=i:i+2*n2
                if ismember(t,sp_fire_time(j,:))
                    last=l4_ext_fire_time(j,:);
                    last=last(last~=0);
                    sp_l4_ext(i,j,t+1)=update_weak(last(end),t,sp_l4_ext(i,j,t),amp_weak,tau_weak);
                else
                    sp_l4_ext(i,j,t+1)=sp_l4_ext(i,j,t);
                end
            end
            for j=1:l4_inh
                if ismember(t,l4_inh_fire_time(j,:))
                    last=l4_inh_fire_time(j,:);
                    last=last(last~=0);
                    inh_ext(i,j,t+1)=update_weak(last(end),t,inh_ext(i,j,t),amp_weak,tau_weak);
                else
                    inh_ext(i,j,t+1)=inh_ext(i,j,t);
                end
            end
            v_l4_ext_input(i,t)=0;
            for j=i:2*n2+i
                v_l4_ext_input(i,t)=v_l4_ext_input(i,t)+sp_l4_ext(i,j)*v_sp(j,t); 
            end
            for j=i:2*n1+i
                v_l4_ext_input(i,t)=v_l4_ext_input(i,t)+ thal_l4_ext(i,j)*g_thal(j+n2,t);
            end
            for j=1:l4_inh
                v_l4_ext_input(i,t)=v_l4_ext_input(i,t)-inh_ext(i,j)*v_l4_inh(i,t);
            end
            if t==1
                if v_l4_ext_input(i,t)>=v_l4_th
                   l4_ext_fire_time(i,k_l4_ext(i,1))=t;
                   v_l4_ext(i,t)=v_spike;
                   for j=i+n2:i+n2+2*n1
                       last=thal_fire_time(j,:);
                       last=last(last~=0);
                       last=last(last<=t);
                       thal_l4_ext(i,j,t+1)=update_strength(last(end),t,thal_l4_ext(i,j,t),amp_strength,tau_strength);
                   end
                   for j=i:i+2*n2
                       last=sp_fire_time(j,:);
                       last=last(last~=0);
                       last=last(last<=t);
                       sp_l4_ext(i,j,t+1)=update_strength(last(end),t,sp_l4_ext(i,j,t),amp_strength,tau_strength);
                   end
                   for j=1:l4_inh
                       last=l4_inh_fire_time(j,:);
                       last=last(last~=0);
                       last=last(last<=t);
                       inh_ext(i,j,t+1)=update_strength(last(end),t,int_ext(i,j,t),amp_strength,tau_strength);
                   end
                   k_l4_ext(i,1)=k_l4_ext(i,1)+1;
                else
                    v_l4_ext(i,t)=v_l4_ext_input(i,t);
                end
            else
                if t-1==l4_ext_fire_time(i,k_l4_ext(i,1))
                    v_l4_ext(i,t)=0;
                else
                    if v_l4_ext_input(i,t)>=v_l4_th
                        l4_ext_fire_time(i,k_l4_ext(i,1))=t;
                        v_l4_ext(i,t)=v_spike;
                        for j=i+n2:i+n2+2*n1
                            last=thal_fire_time(j,:);
                            last=last(last~=0);
                            last=last(last<=t);
                            thal_l4_ext(i,j,t+1)=update_strength(last(end),t,thal_l4_ext(i,j,t),amp_strength,tau_strength);
                        end
                        for j=i:i+2*n2
                            last=sp_fire_time(j,:);
                            last=last(last~=0);
                            last=last(last<=t);
                            sp_l4_ext(i,j,t+1)=update_strength(last(end),t,sp_l4_ext(i,j,t),amp_strength,tau_strength);
                        end
                        for j=1:l4_inh
                            last=l4_inh_fire_time(j,:);
                            last=last(last~=0);
                            last=last(last<=t);
                            inh_ext(i,j,t+1)=update_strength(last(end),t,int_ext(i,j,t),amp_strength,tau_strength);
                        end
                        k_l4_ext(i,1)=k_l4_ext(i,1)+1;
                    else
                        v_l4_ext(i,t)=v_l4_ext_input(i,t)-0.1*v_l4_ext_input(i,t-1);
                    end
                end
            end
        end
    end
    wt_thal_l4_inh(1,time_step).wt=thal_l4_inh(:,:,step/dt);
    wt_thal_l4_ext(1,time_step).wt=thal_l4_ext(:,:,step/dt);
    wt_sp_l4_inh(1,time_step).wt=sp_l4_inh(:,:,step/dt);
    wt_sp_l4_ext(1,time_step).wt=sp_l4_ext(:,:,step/dt);
    wt_inh_ext(1,time_step).wt=inh_ext(:,:,step/dt);
    for i=1:thal
        a=thal_fire_time(i,:);
        a=a(a~=0);
        last_thal(i,1)=a(end)-step/dt;
    end
    for i=1:sp
        a=sp_fire_time(i,:);
        a=a(a~=0);
        last_sp(i,1)=a(end)-step/dt;
    end
    for i=1:l4_inh
        a=l4_inh_fire_time(i,:);
        a=a(a~=0);
        last_l4_inh(i,1)=a(end)-step/dt;
    end
    for i=1:l4_ext
        a=l4_ext_fire_time(i,:);
        a=a(a~=0);
        last_l4_ext(i,1)=a(end)-step/dt;
    end 
    cd('C:\Users\àdmin\Desktop\adarsh\MATLAB\model\running');
    
 end
            
            
                
            
                
            
                    
            
            
            
            
                
                
               
                
                    
            
                
        
                
                
        
