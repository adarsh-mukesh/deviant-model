d1='J:\adarsh_model';
cd('J:\adarsh_model\results\csi_trial_no_std');
ff=dir;
ff(1:2)=[];
n={ff.name};
indx=cellfun(@(v)str2num(v(12:end)),n);
[a,ind]=sort(indx);
ff=ff(ind,:);
bin_sz=0.01;
tot_time=5;
for hh=1%[1 4]
    cd(ff(hh).name)
    fl=dir;fl(1:2)=[];
    load(fl(end).name)
    clearvars -except spont_win iti token_len step fl t_sim d1 step sp_all l4_all ff hh bin_sz tot_time
    fl(end)=[];
    n={fl.name};
    indx=cellfun(@(v)str2num(v(29:end-4)),n);
    [a,ind]=sort(indx);
    fl=fl(ind,:);
    n_toks=t_sim/step;
    figure
    for jj=1:2
        l4_spk_time=cell(100,1);
        sp_spk_time=cell(100,1);
        for ii=1:n_toks
            load(fl((length(fl)/2)*(jj-1)+ii).name);
            l4_spk_time{ii,1}=(l4_fire_time(l4_fire_time>0))./1000;
            sp_spk_time{ii,1}=(sp_fire_time(sp_fire_time>0))./1000;
        end
        xx=[bin_sz*1000:bin_sz*1000:tot_time*1000+bin_sz*1000]./1000;
        subplot(2,1,jj)
        for mm=1:14
            rectangle('Position',[0.5+(mm-1)*0.3 18 0.05 1],'FaceColor','g','EdgeColor','g')
            hold on
        end
        rectangle('Position',[0.5+6*0.3 18 0.05 1],'FaceColor','r','EdgeColor','r')
        hold on
        cd(d1);
        aa=cell2psth(sp_spk_time,bin_sz,tot_time);
        plot(xx,aa,'r')
        hold on
        bb=cell2psth(l4_spk_time,bin_sz,tot_time);
        plot(xx,bb,'b');
        ylim([0 20])
        cd('J:\adarsh_model\results\csi_trial_no_std');
        cd(ff(hh).name)
    end
    cd('J:\adarsh_model\results\csi_trial_no_std');
end
            