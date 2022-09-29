tau_syn=10;%milliseconds (to be put in same units)
t_delay=1;%milliseconds (to be put in same units)
tSim=10;%seconds (to be put in same units)
dt=0.001;%seconds (to be put in same units)
fr=10;%Hz (to be put in same units
n1=8;


for i=1:17
    [spk_mat] = spk_gen_poss(fr,tSim,dt);
    k=1;
    for j=1:length(spk_mat)
        thal_fire_time(i,k)=spk_mat(j);
        k=k+1;
    end
end
for i=1:17
    [g_thal(i,:)]=stim(tau_syn,t_delay,tSim,dt,thal_fire_time(i,:));%thalamic conductances
end


% for i=1:17
%     figure
%     plot(g_thal(i,:));
% end
    weight_init=0.5*exp(-(([1:2*n1+1]-(n1+1)).^2)/(2*((2*n1+2)/4)^2));
    
   input=weight_init*g_thal; 
   
    figure
    plot(input);
    hold on
    plot(g_thal(9,:),'r');
    grid on
    