tic
cd_file=pwd;
result_file='J:\adarsh_model\results';
iti=250;
token_len=50;
dt=0.001;
fs=0.5;
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
step=60*(iti+token_len)/1000;% seconds
u_th=1;
u_sp=1;
amp_strength=0.015;
%amp_strength=0;
amp_weak=0.021;
%amp_weak=0;l4_last_spk
tau_strength=13;%in millisecond
tau_weak=20;%in millisecond
beta=5;
ref_length=20;
t_sim=500*step;
v_th_sp=0.05;
v_th_l4=0.1;
x_thal=zeros(2,3,step*(1/dt)); x_thal(:,1,1)=1;
x_sp=zeros(step*(1/dt),3); x_sp(1,1)=1;
thal_fire_time=zeros(2,500);
sp_fire_time=zeros(1,500);
l4_fire_time=zeros(1,500);
wt_thal_l4=zeros(2,t_sim/step);
wt_sp_l4=zeros(1,t_sim/step);
stop_l4=0;
%init weights
thal_sp=0.2*ones(2,1);
thal_l4=zeros(2,step*(1/dt)); thal_l4(:,1)=0.02;
sp_l4=zeros(1,step*(1/dt)); sp_l4(1,1)=0.22;
iter_num=1;
lps=0.1;
ts=0;
Tir=4200;
tir=Tir;
Tei_thal=10;
Tei_sp=27;
plast=zeros(2,1);
for step_num=1:step:t_sim
    fr=binary_fire(step,token_len,iti,dt);
    thal_fire_time=zeros(2,500);
    if iter_num~=1
        x_thal=zeros(2,3,step*(1/dt)); x_thal(:,:,1)=last_x_thal;
        x_sp=zeros(step*(1/dt),3); x_sp(1,:)=last_x_sp;
        thal_l4=zeros(2,step*(1/dt)); thal_l4(:,1)=last_thal_l4;
        sp_l4=zeros(1,step*(1/dt)); sp_l4(1,1)=last_sp_l4;
    end
    tin=step_num;
    tout=step_num+step;
    for i=1:2
        [spk_mat]=spk_gen_poss_nonhom(fr(i,:),tin,tout,dt);
        thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);
    end
    if step_num~=1;
        thal_fire_time=[thal_fire_time_las, thal_fire_time]; % adding spikes from previous iter
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
    if iter_num~=1
        k_sp=2;
        k_l4=2;
        sp_fire_time(1)=sp_fire_time_las;
        l4_fire_time(:,1)=l4_fire_time_las;
    end
    for t=2:step*(1/dt)
        for i=1:2
            if ismember(t,thal_fire_time(i,:))
                del=1;
            else
                del=0;
            end
            x_thal(i,:,t)=depressing_synapses(x_thal(i,:,t-1),0.9,Tei_thal,tir,del);
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
        x_sp(t,:)=depressing_synapses(x_sp(t-1,:),0.9,Tei_sp,tir,del);
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
        if stop_l4==0
            for i=1:2
                if ismember(t,thal_fire_time(i,:))
                    if l4_last_spk~=0
                        thal_l4(i,t)=update_weak(l4_last_spk,t,thal_l4(i,t-1),amp_weak,tau_weak);
                    else
                        thal_l4(i,t)=thal_l4(i,t-1);
                    end
                else
                    thal_l4(i,t)=thal_l4(i,t-1);
                end
            end
        else
            for i=1:2
                thal_l4(i,t)=thal_l4(i,t-1);
            end
        end
        if ismember(t,sp_fire_time)
            if 1%sum(plast)==2 
                sp_l4(1,t)=update_weak(l4_last_spk,t,sp_l4(1,t-1),amp_weak,tau_weak);
            else
                sp_l4(1,t)=sp_l4(1,t-1);
            end
        else
            sp_l4(1,t)=sp_l4(1,t-1);
        end
        if l4_last_spk==0
            v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,t-1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,t-1)*g_sp(1,t-1));
        else
            if t-l4_last_spk<=ref_length
                v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,t-1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,t-1)*g_sp(1,t-1))-beta*exp(-(t-l4_last_spk-1)/t_refr);
            else
                v_l4(1,t)=(1-lps)*((squeeze(x_thal(:,2,t-1)).*thal_l4(:,t-1))'*g_thal(:,t-1)+x_sp(t-1,2)*sp_l4(1,t-1)*g_sp(1,t-1));
            end
        end
        if ((v_l4(1,t)>v_th_l4) || (rand<=1-exp(-1*fs*dt)))
            if stop_l4==0
                for i=1:2
                    up_lim=0.3;
                    thal_l4(i,t)=update_strength(thal_fire_time(i,:),t,thal_l4(i,t-1),amp_strength,amp_weak,tau_strength,up_lim);
                    if thal_l4(i,t)==up_lim
                        stop_l4=0;
                        plast(i)=1;
                    end
                end
            else
                for i=1:2
                    thal_l4(i,t)=thal_l4(i,t-1);
                end
            end
            up_lim=0.3;
            sp_l4(1,t)=update_strength(sp_fire_time,t,sp_l4(1,t-1),amp_strength,amp_weak,tau_strength,up_lim);
            if sum(plast)~=2
                if sp_l4(1,t)<sp_l4(1,t-1)
                    sp_l4(1,t)=sp_l4(1,t-1);
                end
            end
%             if sum(plast)==2 
%                 sp_l4(1,t)=update_strength(sp_fire_time,t,sp_l4(1,t-1),amp_strength,amp_weak,tau_strength,up_lim);
%             else
%                 sp_l4(1,t)=sp_l4(1,t-1);
%             end
            %sp_l4(1,t)=update_strength(sp_fire_time,t,sp_l4(1,t-1),amp_strength,amp_weak,tau_strength,up_lim);
            %sp_l4(1,t)=sp_l4(1,t-1);
            l4_last_spk=t;
            l4_fire_time(1,k_l4)=t;
            k_l4=k_l4+1;
            del=1;
        end
    end
    thal_fire_time_las=zeros(2,1);
    for i=1:2
        a=thal_fire_time(i,:);
        a=a(a>0);
        if isempty(a)~=1
            thal_fire_time_las(i,1)=a(end)-step*(1/dt);
        end
    end
    sp_fire_time_las=sp_last_spk-step*(1/dt);
    l4_fire_time_las=l4_last_spk-step*(1/dt);
    last_x_thal=x_thal(:,:,step*(1/dt));
    last_x_sp=x_sp(step*(1/dt),:);
    last_thal_l4=thal_l4(:,end);
    wt_thal_l4(:,iter_num)=thal_l4(:,end);
    last_sp_l4=sp_l4(1,end);
    wt_sp_l4(1,iter_num)=sp_l4(1,end);
    cd(result_file);
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(iter_num));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_fire_time sp_l4 thal_l4 g_thal g_sp v_l4 v_sp x_thal x_sp fr tir',resultfile));
    cd(cd_file);
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
%    if sum(plast)>0
%        tir=tir*0.997;
%    end
   
end
cd(result_file);
save sum
cd(cd_file)
figure
plot(wt_thal_l4(2,:));
hold on
plot(wt_thal_l4(1,:),'r');
hold on
plot(wt_sp_l4(1,:),'k');

toc