function [fr]=fire_ssa_gen(f1,f2,token_len,iti,num_token,dev_pos,tun_wid,gap)
global sp n1
if dev_pos>num_token
    error('devint position must be less then token number');
end
a=10*exp(-(([1:25]-f1).^2)/(2*tun_wid^2));
fr=repmat([repmat(max(a,1)',1,token_len) ones(25,iti)],1,15);
fr=[fr ones(25,gap)];
fr(:,(dev_pos-1)*(token_len+iti)+1:(dev_pos-1)*(token_len+iti)+token_len)=repmat((max(10*exp(-(([1:25]-f2).^2)/(2*tun_wid^2)),1))',1,token_len);
end


        