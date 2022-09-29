function [new_weight]=update_weak(last_fire,t,old_weight,amp_weak,tau_weak)

del_t=t-last_fire;
new_weight=old_weight*(1-amp_weak*exp(-1*abs(del_t)/tau_weak));
end
