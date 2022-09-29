function gen_psth=cell2psth(dat,bin_sz,total_time)
%%% bin size in seconds
%%% all the iters should be in form of a verticle cell 
%%% size(dat,1)=number of iters; size(dat,2)=1; 
nreps=size(dat,1);
nbins=fix(total_time/bin_sz)+1;
psth=zeros(nreps,nbins);
for ii=1:size(dat,1)
    spks=dat{ii,1};
    nspks=length(spks);
    for jj=1:nspks
        psth(ii,fix(spks(jj)/bin_sz)+1)=psth(ii,fix(spks(jj)/bin_sz)+1)+1;
    end
end
gen_psth=mean(psth)/(bin_sz);
end