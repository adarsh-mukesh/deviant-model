function[alpha]=waiting(evn_t)
evn_t=evn_t(evn_t~=0);
k=length(evn_t);
for i=1:k-1
    alpha(i)=evn_t(i+1)-evn_t(i);
end
figure
hist(alpha);
    

end