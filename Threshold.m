function [spikeMatrix, tVec] = Threshold(population,value)
for j=1:population
    [spikeMat,tVec]=poissonSpikeGen(30,1,20);
    spikeMat_sum(:,:,j)=spikeMat;
end
for k=1:20
    for l=1:1000
        if sum(spikeMat_sum(k,l,:))>=value
            spikeMatrix(k,l)=1;
        else
            spikeMatrix(k,l)=0;
        end
    end
end
end






