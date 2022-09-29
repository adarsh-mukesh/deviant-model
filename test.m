for t=2:18000
    a=sp_fire_time(sp_fire_time>t-(10*tau_syn));
    a=a(a<t);
    a=a(a>0);
    if isempty(a)
        g(t)=0;
    else
        g(t)=0;
        for nn=1:length(a)
            g(t)=g(t)+exp(-(t-a(nn))/tau_syn);
        end
    end
end