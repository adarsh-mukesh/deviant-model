
% declare parameters
dt=1;
t_end=1000;
t_stimstart=200;
t_stimend=800;
E_L=-70;
V_th=-63;
V_reset=-75;
V_spike=20;
R_m=10;
tau=10;

%initial value
t_vect=1:dt:t_end;
V_vect=zeros(1,length(t_vect));
i=1;
V_vect(i)=E_L;
V_plot_vect(i)=V_vect(i);

%non-hom inj
%I_e_vect=k;
%for i=1:199
 %   I_e_vect(i)=0;
%end
%for i=801:1000
 %  I_e_vect(i)=0;
%end
%i=1;
%current stimulus
I_stim=1;
I_e_vect=zeros(1,t_stimstart/dt);
I_e_vect=[I_e_vect I_stim*ones(1,1+((t_stimend-t_stimstart)/dt))];
I_e_vect=[I_e_vect zeros(1,(t_end-t_stimend)/dt)];
numspikes=0;

%integrate and fire

for t=dt:dt:t_end-1
     V_inf=E_L+I_e_vect(i)*R_m;
     V_vect(i+1)=V_inf+(V_vect(i)-V_inf)*exp(-dt/tau);
     if V_vect(i+1)>V_th
         V_vect(i+1)=V_reset;
         V_plot_vect(i+1)=V_spike;
         numspikes=numspikes+1;
     else
         V_plot_vect(i+1)=V_vect(i+1);
     end
     i=i+1;
end
avgrate=1000*numspikes/(t_stimend-t_stimstart);

%plots
figure(1)
plot(t_vect,V_vect);
title('Voltage vs. Time');
xlabel('Time in ms');
ylabel('Voltage in mV');

figure(2)
plot(t_vect,V_plot_vect);
title('Spike vs. Time');
xlabel('Time in ms');
ylabel('Voltage in mV');
hold on

