% to simulate depressing synapses, xrecovery xeffective xinactive are
% neurotransmitter concentration of a neuron that changes dynamically with time.
%Tre, Tei, Tir are the respective time constants. 
% delta=1 for spiking event and delta=0 for no spike 
function [x_new]=depressing_synapses(x_old,Tre,Tei,Tir,delta)
x_recovery=x_old(1,1,1);
x_effective=x_old(1,2,1);
x_inactive=x_old(1,3,1);

ratex_recovery=-delta*x_recovery/Tre + x_inactive/Tir;
ratex_effective=delta*x_recovery/Tre - x_effective/Tei;

x_recovery=x_recovery+ratex_recovery;
x_effective=x_effective+ratex_effective;    
    if x_recovery<0
        x_recovery=0;
        
    elseif x_recovery>1
        x_recovery=1;
        
    end
    
    if x_effective<0
        x_effective=0;
    
    elseif x_effective>1;
        x_effective=1;
        
    end
    
    x_inactive=1-(x_effective+x_recovery);
    if x_inactive<0
        x_inactive=0;
    elseif x_inactive>1
        x_inactive=1;
    end
    x_new=[x_recovery x_effective x_inactive];
        
end 