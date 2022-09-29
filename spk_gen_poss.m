function [spk_mat] = spk_gen_poss(fr,tin,tout,dt)
spk_mat=0;
tSim=tout-tin;
nBins=floor(tSim/dt);
rng('shuffle');
for i=1:nBins
    if rand<=1-exp(-1*fr*dt)
        spikeMat(i)=1;
    else
        spikeMat(i)=0;
    end
end
tVec=0:dt:tSim-dt;
k=1;
for i=1:size(spikeMat,2)
    if spikeMat(i)==1;
        spk_mat(k)=i;
        k=k+1;
    end
end

        
        




end

