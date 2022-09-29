function [fr]=fire_oddball_seq(f1,f2,token_len,iti,tun_wid,step,dt)
% step shpuld be an integral multiple of (iti+token_len)
global sp n1
fr=zeros(2*n1+sp,step/dt);
for t=1:(iti+token_len):step/dt
    rng('shuffle');
    if rand>0.9
        fr(:,t:t+iti+token_len-1)=[repmat((max(8*exp(-(([1:2*n1+sp]-f2).^2)/(2*tun_wid^2)),0.1))',1,token_len) 0.1*ones(2*n1+sp,iti)];
    else
        fr(:,t:t+iti+token_len-1)=[repmat((max(8*exp(-(([1:2*n1+sp]-f1).^2)/(2*tun_wid^2)),0.1))',1,token_len) 0.1*ones(2*n1+sp,iti)];
    end
end
