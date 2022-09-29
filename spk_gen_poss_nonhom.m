function [spk_mat] = spk_gen_poss_nonhom(fr,tin,tout,dt)
spk_mat=[];
tSim=tout-tin;
nBins=floor(tSim/dt);
rng('shuffle');
for i=1:nBins
    if rand<=1-exp(-1*fr(1,i)*dt);
        spikeMat(i)=1;
    else
        spikeMat(i)=0;
    end
end
spk_mat=find(spikeMat==1);
end

