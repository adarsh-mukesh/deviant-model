function [fr]=fire_allfreq_seq(token_len,iti,tun_wid,step,dt)
%step shpuld be an integral multiple of (iti+token_len)
global sp n1
%sp=33; n1=8;
fr=zeros(2*n1+sp,step/dt);
for t=1:(iti+token_len):step/dt
    rng('shuffle');
    fr(:,t:t+iti+token_len-1)=[repmat((max(10*exp(-(([1:2*n1+sp]-randi([1 41],1,1)).^2)/(2*tun_wid^2)),0.1))',1,token_len) 0.1*ones(2*n1+sp,iti)];
end