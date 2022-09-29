function [new_weight] = update_strength_inhib(fire_time,t,old_weight,amp_strength,amp_weak,tau_strength,up_lim)
fire_time=fire_time(fire_time>0);
c=fire_time(fire_time<=t);
if isempty(c)==1
    new_weight=old_weight;
else
    del_t=t-c(end);
    new_weight=old_weight*(1+amp_strength*exp(-1*abs(del_t)/tau_strength));
    if del_t==0;
        new_weight=old_weight*(1-amp_weak);
    end
    if new_weight>up_lim
        new_weight=up_lim;
    end
end

end

