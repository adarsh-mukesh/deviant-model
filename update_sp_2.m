function[new_weight]=update_sp_2(old_weight,t,f_t)
del_t=f_t-t;
new_weight=old_weight*(1-(0.0001*exp((-abs(del_t))/20)));
end