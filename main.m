tic
% defining parameters-----------------------------
dt=0.001;
tSim=10;
fr=20;
sp=48;
n1=5;
n2=8;
alpha_tau=0.02;
E_L=-70;
%------------------------------------------------

% weight generation----------------------------------
l4=sp-2*n2;
thal_sp=zeros(sp,sp+2*n1,tSim/dt);
thal_l4=zeros(l4,l4+2*n1,tSim/dt);
sp_l4=zeros(l4,sp,tSim/dt);
[thal_sp(:,:,1), thal_l4(:,:,1), sp_l4(:,:,1)]=weight_gen(sp,n1,n2);
%----------------------------------------------------

%stim generation--------------------------------------
for i=1:sp+2*n1
    [a]=spk_gen_poss(fr,tSim,dt);
    k=1;
    for j=1:size(a,2)
        fire_time(i,k)=a(j);
        k=k+1;
    end
end
lim=size(fire_time,2);
for i=1:sp+2*n1
    [b]=stim(alpha_tau,tSim,dt,fire_time(i,:));
    k=1;
    for j=1:size(b,2)
        inject(i,k)=b(j);
        k=k+1;
    end
    for j=size(b,2)+1:tSim/dt
        inject(i,j)=0;
    end
end
%-----------------------------------------------------
l4=sp-2*n2;
% index init------------------------------------------
k_sp=ones(sp,1);
k_l4=ones(sp-2*n2,1);
sp_fire_las=zeros(sp,1);
l4_fire_las=zeros(sp,1);
%------------------------------------------------------

% voltage generation-----------------------------------
sp_volt(1:sp,1)=E_L;
l4_volt(1:sp-2*n2,1)=E_L;
%------------------------------------------------------
count1=0;
%run thal to sp--------------------------------------------------
for t=1:tSim/dt-1
    for i=1:sp
% plasticity thal firing-------------------------------------------------
        if sp_fire_las(i,1)~=0
            for j=i:i+2*n1
                k=1;
                while k<lim
                    if fire_time(j,k)==0
                        break
                    end
                    if fire_time(j,k)==t
                        count2=1;
                        [thal_sp(i,j,t+1)]=update_thal(thal_sp(i,j,t),t,sp_fire_las(i,1));
                    end
                    k=k+1;
                end
            end
        end
%-------------------------------------------------------------------------
        sp_stim=0;
        for j=1:2*n1+1
            sp_stim=sp_stim+inject(i+j-1,t)*thal_sp(i,i+j-1,t);
        end
        [a,b]=sb_int_n_fire(sp_stim,sp_volt(i,t));
        sp_volt(i,t+1)=b;
        if a==1
            count3=5;
            sp_fire(i,k_sp(i))=t;
            sp_fire_las(i,1)=t;
            k_sp(i)=k_sp(i)+1;
% plasticity sp fire----------------------------------------------------
            for j=i:i+2*n1
                count4=6;
                thal_sp(i,j,t+1)=update_sp(thal_sp(i,j,t),t,fire_time(j,:));
                if t+1~=1;
                    count4=99;
                end
            end
%-----------------------------------------------------------------------
           
        end
        if a==0
            for j=i:i+2*n1
                cn=0;
                k=1;
                while k<lim
                    if fire_time(j,k)==t
                        if sp_fire_las(i,1)==0
                            thal_sp(i,j,t+1)=thal_sp(i,j,t);
                        end
                        cn=cn+1;
                    end
                    k=k+1;
                end
                if cn==0;
                    thal_sp(i,j,t+1)=thal_sp(i,j,t);
                end
            end
        end
    end
end
%sp to l4---------------------------------------------------------
for ii=1:sp
    [b]=stim(alpha_tau,tSim,dt,sp_fire(ii,:));
    k=1;
    for jj=1:size(b,2)
        inject_sp(ii,k)=b(jj);
        k=k+1;
    end
end
  
%------------------------------------------------------------------
lim1=size(sp_fire,2);
% current injection to l4------------------------------------------
for t=1:tSim/dt-1
    for i=1:l4
% plasticity thal and sp firing------------------------------------        
        if l4_fire_las(i,1)~=0
            for j=i:i+2*n2
                k=1;
                while k<lim1
                    if sp_fire(j,k)==0
                        break
                    end
                    if sp_fire(j,k)==t
                        sp_l4(i,j,t+1)=update_sp_2(sp_l4(i,j,t),t,l4_fire_las(i,1));
                    end
                    k=k+1;
                end
            end
            for j=i+n2:i+n2+2*n1
                k=1;
                while k<lim
                    if fire_time(j,k)==0
                        break
                    end
                    if fire_time(j,k)==t
                        thal_l4(i,j-n2,t+1)=update_thal(thal_l4(i,j-n2,t),t,l4_fire_las(i,1));
                    end
                    k=k+1;
                end
            end
        end                
%----------------------------------------------------------------------------------        
        l4_stim=0;
        for ii=1:2*n1+1
            l4_stim=l4_stim+inject(i+ii-1+n2,t)*thal_l4(i,i+ii-1);
        end
        for ii=1:2*n2+1
            l4_stim=l4_stim+inject_sp(i+ii-1,t)*sp_l4(i,i+ii-1);
        end
        [a,b]=sb_int_n_fire(l4_stim,l4_volt(i,t));
        l4_volt(i,t+1)=b;
        if a==1
            l4_fire(i,k_l4(i))=t;
            l4_fire_las(i,1)=t;
            k_l4(i)=k_l4(i)+1;
% plasticity l4 firing-------------------------------------------------            
            for j=i:i+2*n2
                sp_l4(i,j,t+1)=update_l4(sp_l4(i,j,t),t,sp_fire(j,:));
            end
            for j=i+n2:i+n2+2*n1
                thal_l4(i,j-n2,t+1)=update_l4(thal_l4(i,j-n2,t),t,fire_time(j,:));
            end
%--------------------------------------------------------------------------            
        end
        if a==0
            for j=i:i+2*n2
                k=1;
                cn=0;
                while k<lim1
                    if sp_fire(j,k)==t
                        if l4_fire_las(i,1)==0
                            sp_l4(i,j,t+1)=sp_l4(i,j,t);
                        end
                        cn=cn+1;
                    end
                    k=k+1;
                end
                if cn==0
                    sp_l4(i,j,t+1)=sp_l4(i,j,t);
                end
            end
            for j=i+n2:i+n2+2*n1
                k=1;
                cn=0;
                while k<lim
                    if fire_time(j,k)==t
                        if l4_fire_las(i,1)==0
                            thal_l4(i,j-n2,t+1)=thal_l4(i,j-n2,t);
                        end
                        cn=cn+1;
                    end
                    k=k+1;
                end
                if cn==0
                    thal_l4(i,j-n2,t+1)=thal_l4(i,j-n2,t);
                end
            end
        end         
    end
end
toc
            
          
    




        