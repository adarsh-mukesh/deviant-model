function [firing,V_next]=sb_int_n_fire_sp(curr,V)
 E_L=-70;
 R_m=10;
 dt=1;
 tau=10;
 V_th=-55;
 V_reset=-75;
 if V==V_reset
     curr=0;
 end  
 dV=(dt/tau)*(R_m*curr-(V-E_L));
 V_next=V+dV;
 if V_next>V_th
     firing=1;
     V_next=V_reset;
 else 
     firing=0;
 end
 end