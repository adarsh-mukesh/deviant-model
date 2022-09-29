tic
%declare_parameters--------------------------------
global sp n1 n2 n3
sp=17;
n1=4;
n2=4;
n3=4;
sigma=1;
up_lim=0.2; 
tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
t_refr=2;%milliseconds (to be put in same units)
% stim params-------------------------------
f1=13;
f2=10;
token_len=50;%milliseconds
iti=250;%milliseconds
num_token=15;
dev_pos=9;
tun_wid=1;
gap=000;%milliseconds
%-------------------------------------------
%step=10;
step=(15*(token_len+iti)+gap)/1000;%time in one step (in seconds)
t_Sim=100*step;% in seconds ; total simulation time seconds (to be put in same units)
dt=0.001;%seconds (to be put in same units)
fr=fire_ssa_gen(f1,f2,token_len,iti,num_token,dev_pos,tun_wid,gap);%Hz (to be put in same units)
fs=0.5;
%fs_w=0.2;
lps=0.1;% leak par sample
v_th_sp=0.03;
v_th_l4_inh=0.07;
v_th_l4_ex=0.07;
u_th=1;
u_sp=1;
u_l4=1;
v_spike=10;
%amp_strength=0.015;
amp_strength=0;
%amp_weak=0.021;
amp_weak=0;
tau_strength=13;%in millisecond
tau_weak=20;%in millisecond
beta=5;
ref_length=20;
Tre=0.9;%in milliseconds
Tei=5.3;%in milliseconds
Tir=5000;%in milliseconds
%--------------------------------------------------
%weight_generation---------------------------------
l4_inh=sp-2*n2;
l4_ex=l4_inh-2*n3;
thal=sp+2*n1;
thal_l4_inh=zeros(l4_inh,l4_inh+2*n1);
thal_l4_ex=zeros(l4_ex,l4_ex+2*n1);
sp_l4_inh=zeros(l4_inh,sp);
sp_l4_ex=zeros(l4_ex,l4_ex+2*n2);
%l4_inh_ex=zeros(l4_ex,l4_inh,step/dt);
[thal_sp, thal_l4_inh, sp_l4_inh, thal_l4_ex, sp_l4_ex, l4_inh_ex]=weight_gen_abcd(sigma);% weight matrix
%thal_l4(:,:,1)=a;
%--------------------------------------------------------------
thal_fire_time=[];
sp_fire_time=zeros(sp,1);
l4_inh_fire_time=zeros(l4_inh,1);
l4_ex_fire_time=zeros(l4_ex,1);
iter_num=1;
for step_num=1:step:t_Sim
    alp=[];
    %fr=50*exp(-(([1:25]-iter_num).^2)/(2*tun_wid^2));
    x_thal=zeros(thal,3,step*(1/dt)); x_thal(:,1,1)=1;
    x_sp=zeros(sp,3,step*(1/dt)); x_sp(:,1,1)=1;
    x_l4_inh=zeros(l4_inh,3,step*(1/dt)); x_l4_inh(:,1,1)=1;
    %stim generation-----------------------------------------------
    tin=step_num;
    tout=step_num+step;
    g_thal=zeros(thal,step/dt);
    alp=[];
    thal_fire_time=zeros(thal,500);
    for i=1:thal
        [spk_mat] = spk_gen_poss_nonhom(fr(i,:),tin,tout,dt);
        alp=[alp;fr(i,:)];
        %[spk_mat]=spk_gen_poss(fr(i),tin,tout,dt);
        thal_fire_time(i,1:length(spk_mat))=sort(spk_mat);
        [g_thal(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,thal_fire_time(i,:),u_th);%thalamic conductances
    end
    %g_thal=ones(thal,step/dt);
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
            x_thal(i,:,t)=depressing_synapses(x_thal(i,:,t-1),0.9,5.3,2000,del);
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
            x_sp(i,:,t)=depressing_synapses(x_sp(i,:,t-1),0.9,20,2000,del);
        end
    end
%         for i=1:sp
%             [spont]=spk_gen_poss(fs,tin,tout,dt);
%             a=sp_fire_time(i,:);
%             a=a(a>0);
%             a=sort([a,spont]);
%             sp_fire_time(i,:)=0;
%             for j=1:length(a)
%                 sp_fire_time(i,j)=a(j);
%             end
%         end
%         a=[];
    fprintf('thal and sp firing complete\n');
    for i=1:sp
        [g_sp(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,sp_fire_time(i,:),u_sp);
    end
    %g_sp=ones(sp,step/dt);
    sp_fire_time(:,~any(sp_fire_time,1))=[];
    %-------------------------------------------------------------------
    k=ones(l4_inh,1);
    g_l4_inh=zeros(l4_inh,step/dt);
    l4_inh_last_spk=zeros(l4_inh,1);% index of the l4_fire_time m  atrix for each l4 neuron
    v_l4_inh=zeros(l4_inh,step/dt);
    l4_inh_fire_time=zeros(l4_inh,500);
    for t=2:step/dt
        for i=1:l4_inh
            % now check which sp and thal neurons have fired on this time
            % and see the difference from the last time when the
            % corresponding l4 fired
            
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
                del=1;
            else
                del=0;
            end
            %x_l4_inh(i,:,t)=depressing_synapses(x_l4_inh(i,:,t-1),Tre,Tei,Tir,del);
            x_l4_inh(i,:,t)=x_l4_inh(i,:,t-1);
        end
    end
%     for i=1:l4_inh
%         [spont]=spk_gen_poss(fs,tin,tout,dt);
%         a=l4_inh_fire_time(i,:);
%         a=a(a>0);
%         a=sort([a,spont]);
%         l4_inh_fire_time(i,:)=0;
%         for j=1:length(a)
%             l4_inh_fire_time(i,j)=a(j);
%         end
%     end
%     a=[];
    
    l4_inh_fire_time(:,~any(l4_inh_fire_time,1))=[];
    for i=1:l4_inh
        [g_l4_inh(i,:)]=stim(tau_syn,t_delay,tin,tout,dt,l4_inh_fire_time(i,:),u_l4);
    end
    k=ones(l4_ex,1);
    l4_ex_last_spk=zeros(l4_ex,1);% index of the l4_fire_time m  atrix for ]4 neuron
    v_l4_ex=zeros(l4_ex,step/dt);
    l4_ex_fire_time=zeros(l4_ex,500);
    for t=2:step/dt
        for i=1:l4_ex
            if l4_ex_last_spk(i,1)==0
                v_l4_ex(i,t)=(1-lps)*((squeeze(x_thal(i+n2+n3:i+2*n1+n2+n3,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+(squeeze(x_sp(i+n3:i+n3+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i+n3:i+n3+2*n2,t-1)-(squeeze(x_l4_inh(i:i+2*n3,2,t-1))'.*l4_inh_ex(i,i:i+2*n3))*g_l4_inh(i:i+2*n3,t-1));
            else
                if t-l4_ex_last_spk(i,1)<=ref_length
                    v_l4_ex(i,t)=(1-lps)*((squeeze(x_thal(i+n2+n3:i+2*n1+n2+n3,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+(squeeze(x_sp(i+n3:i+n3+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i+n3:i+n3+2*n2,t-1)-(squeeze(x_l4_inh(i:i+2*n3,2,t-1))'.*l4_inh_ex(i,i:i+2*n3))*g_l4_inh(i:i+2*n3,t-1))-beta*exp(-(t-l4_ex_last_spk(i,1)-1)/t_refr);
                else
                    v_l4_ex(i,t)=(1-lps)*((squeeze(x_thal(i+n2+n3:i+2*n1+n2+n3,2,t-1))'.*thal_l4_inh(i,i:i+2*n1))*g_thal(i+n2+n3:i+n2+n3+2*n1,t-1)+(squeeze(x_sp(i+n3:i+n3+2*n2,2,t-1))'.*sp_l4_inh(i,i:i+2*n2))*g_sp(i+n3:i+n3+2*n2,t-1)-(squeeze(x_l4_inh(i:i+2*n3,2,t-1))'.*l4_inh_ex(i,i:i+2*n3))*g_l4_inh(i:i+2*n3,t-1));
                end
            end
            if v_l4_ex(i,t)>=v_th_l4_ex
                l4_ex_last_spk(i,1)=t;
                l4_ex_fire_time(i,k(i,1))=t;
                k(i,1)=k(i,1)+1;
            end
        end
    end
    
%     for i=1:l4_ex
%         [spont]=spk_gen_poss(fs,tin,tout,dt);
%         a=l4_ex_fire_time(i,:);
%         a=a(a>0);
%         a=sort([a,spont]);
%         l4_ex_fire_time(i,:)=0;
%         for j=1:length(a)
%             l4_ex_fire_time(i,j)=a(j);
%         end
%     end
%     a=[];
    cd('K:\adarsh_model\results');
    str=date;
    resultfile=strcat('result_',str,'_iter_',num2str(iter_num));
    eval(sprintf('save %s thal_fire_time sp_fire_time l4_inh_fire_time l4_ex_fire_time thal_sp sp_l4_inh sp_l4_ex thal_l4_inh thal_l4_ex g_thal g_sp g_l4_inh v_l4_inh v_l4_ex v_sp x_thal x_sp x_l4_inh',resultfile))
    cd('K:\adarsh_model');
    
    fprintf('iter %i complete\n',iter_num);
    iter_num=iter_num+1;
end
cd('K:\adarsh_model\results');
save sum
cd('K:\adarsh_model');
toc