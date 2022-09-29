d1='J:\adarsh_model\results\all_early\bin_sd_early_f5';
cd(d1);
f1=dir('*.mat');f1(end)=[];
[~,idx]=sort(cellfun(@(w) str2num(w(max(find(w=='_'))+1:end-4)),{f1.name}));
f1=f1(idx);
htt=cell(1,2);
upm=50;
for ii=1:250
    load(f1(ii).name)
    l4_fire_time=l4_fire_time(l4_fire_time>0);
    sp_fire_time=sp_fire_time(sp_fire_time>0);
    for jj=1:2
        ap=thal_fire_time(jj,:);
        ap=ap(ap>0);
        for tt=1:length(ap)
            op=ap(tt)-l4_fire_time;
            op=op(op>=1); op=op(op<upm);
            htt{1,jj}=[htt{1,jj}; -1*min(op)];
        end
        for tt=1:length(l4_fire_time)
            op=l4_fire_time(tt)-ap;
            op=op(op>=0); op=op(op<upm);
            htt{1,jj}=[htt{1,jj}; min(op)];
        end
    end
end
htt_all=cellfun(@(w,x) [[w;NaN*ones(max([0 length(x)-length(w)]),1)] [x;NaN*ones(max([0 length(w)-length(x)]),1)]],htt_early,htt_later,'UniformOutput',false);
bn=40;
bsz=5;
binsp=(1:bsz:upm)'; binsp=[binsp binsp+(bsz-1)]; binsn=[0:-1*bsz:-1*upm]'; binsn=[binsn(1:end-1) binsn(1:end-1)-1*(bsz-1)];
bins=[binsn(end:-1:1,[2 1]);binsp];
bar_all=cell(1,2);
for ii=1:2
    ap=htt_all{1,ii};
    a1=ap(:,1); a2=ap(:,2);
    a1=a1(~isnan(a1)); a2=a2(~isnan(a2));
    for jj=1:size(bins,1)
        bar_all{1,ii}(jj,1)=length(find(a1>=bins(jj,1) & a1<=bins(jj,2)))/length(a1);
        bar_all{1,ii}(jj,2)=length(find(a2>=bins(jj,1) & a2<=bins(jj,2)))/length(a2);
    end
end
pp=[]; hh=[];
sign={'left','right','left'};
figure
for ii=1:2
    %figure
    subplot(1,2,ii)
    bar(bins(:,1),bar_all{1,ii})
    ylim([0 0.9])
    hold on
    plot([nanmean(htt_all{1,ii}(:,1)) nanmean(htt_all{1,ii}(:,1))],[0 0.9],'b');
    hold on
    plot([nanmean(htt_all{1,ii}(:,2)) nanmean(htt_all{1,ii}(:,2))],[0 0.9],'y');
    hold on
    plot([0 0],[0 0.9],'k')
    title(nanmean(htt_all{1,ii}));
    ap=htt_all{1,ii}(:,1);
    ap=ap(~isnan(ap));
    bp=htt_all{1,ii}(:,2);
    bp=bp(~isnan(bp));
    [h p]=ttest2(ap,bp,0.05,sign{ii});
    hh=[hh h]; pp=[pp p];
end

        
            
            