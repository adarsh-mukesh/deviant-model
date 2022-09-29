function [fr]=pure_ssatone(f1,f2,token_len,iti,num_token,dev_pos)
global sp n1
if dev_pos>num_token
    error('devint position must be less then token number');
end
a=zeros(sp+2*n1,token_len);
a(f1,:)=50;
a=[a zeros(sp+2*n1,iti)];
fr=repmat(a,1,15);
fr=[fr zeros(57,2000)];
fr(f1,(iti+token_len)*(dev_pos-1):(iti+token_len)*(dev_pos-1)+token_len)=0;
fr(f2,(iti+token_len)*(dev_pos-1):(token_len+iti)*(dev_pos-1)+token_len)=50;

