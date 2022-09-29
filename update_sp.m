function [new_weight]=update_sp(old_weight,t,fire_time)


k=1;
while k<size(fire_time,2)+1
      if fire_time(1,k)==0;
         t_prev=fire_time(1,k-1);
      end
      if fire_time(1,k)>t
          if k-1==0
              new_weight=old_weight;
              return
          end
         t_prev=fire_time(1,k-1);
         break
      end
      k=k+1;
      if k==size(fire_time,2)
         t_prev=fire_time(1,k);
      end
      
end
 del_t=t-t_prev;
 new_weight=old_weight*(1+0.0001*exp((-1*del_t)/3));
 
    
end
