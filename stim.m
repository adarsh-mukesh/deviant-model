function [g]=stim(tau_syn,t_delay,tin,tout,dt,spk_mat,u)
tSim=tout-tin;
spk_mat=spk_mat(spk_mat~=0);
spk_trn=zeros(1,uint16(tSim*(1/dt)));
spk_trn(spk_mat+1)=1;
tVec=1:10*tau_syn;
ker=[zeros(1,t_delay) exp(-(tVec)/tau_syn)];
g=conv(spk_trn,ker);
if length(g)>uint16(tSim*(1/dt))
    g=g(1:uint16(tSim*(1/dt)));
elseif length(g)<uint16(tSim*(1/dt))
        g(length(g):uint16(tSim*(1/dt)))=0;
end
g=u*g;
end

 