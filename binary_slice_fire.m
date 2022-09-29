function [ fr ] = binary_slice_fire(spont_win,token_len,iti,num_tok,dev_loc)
fr1=0.1*ones(1,spont_win);
fr2=fr1;
fr1=[fr1 repmat([8*ones(1,token_len) 0.1*ones(1,iti)],1,num_tok)];
fr2=[fr2 repmat([0.15*ones(1,token_len) 0.1*ones(1,iti)],1,num_tok)];
fr1(spont_win+(dev_loc-1)*(token_len+iti)+1:spont_win+(dev_loc-1)*(token_len+iti)+50)=0.15;
fr2(spont_win+(dev_loc-1)*(token_len+iti)+1:spont_win+(dev_loc-1)*(token_len+iti)+50)=8;
fr=[fr1;fr2];
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


end

