function [fr]=binary_fire(step,token_len,iti,dt)
fr=0.1*ones(2,step/dt);
rng('shuffle');
for i=1:(iti+token_len):step/dt
    if rand>0.9
        fr(2,i:i+token_len-1)=0.15;
        fr(1,i:i+token_len-1)=8;
    else
        fr(1,i:i+token_len-1)=8;
        fr(2,i:i+token_len-1)=0.15;
    end
end