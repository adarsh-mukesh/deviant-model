% This code euns the network saved and brought to equilibrium by
% TY_Loebel_Par_All_Prots_Nevo, according to the specified nev_cond
tic
n_stim    = 100; % Total no. of stimuli (Best take a product of 10)
if Probs_F2L < 0.1
    n_stim = 10/Probs_F2L;
end
t_prot    = n_stim*(duration + ISI) + 2*post_stim;
tmax_tot  = t_eq + t_prot; % Maximum time that the simulation will reach (in seconds) 
num_steps = floor(tmax_tot/dt); % Total number of steps in the simulation

Freq_Scale = 1:M;

Prob.L       = Probs_F2L;
Freqs_Pres.L = [Rec_Column - Freq_Diff/2; Rec_Column + Freq_Diff/2]; 
    
Prob.H       = 1 - Probs_F2L;
Freqs_Pres.H = [Rec_Column + Freq_Diff/2; Rec_Column - Freq_Diff/2]; 

Prob.E       = 0.5;
Freqs_Pres.E = [Rec_Column - Freq_Diff/2; Rec_Column + Freq_Diff/2];
        
Prob.DB = Probs_F2L;
if Prob.DB == 0.1
    flank_freq    = (1/Probs_F2L - 2)/2;
    Freqs_Pres.DB = sort([(Rec_Column - Freq_Diff/2) - (0:Freq_Diff:flank_freq*Freq_Diff), (Rec_Column + Freq_Diff/2) + (0:Freq_Diff:flank_freq*Freq_Diff)]);
end   
    
Prob.DN = Probs_F2L;
if Prob.DN == 0.1
    freq_res      = Freq_Diff*2*Probs_F2L;
    Freqs_Pres.DN = (Rec_Column - Freq_Diff/2) + (-2*freq_res:freq_res:7*freq_res);
end

Prob.DA1 = Probs_F2L;
Freqs_Pres.DA1 = Rec_Column - Freq_Diff/2;
        
Prob.DA2 = Probs_F2L;
Freqs_Pres.DA2 = Rec_Column + Freq_Diff/2;

if Prob.DN == 0.1
    M_aug      = M + sum(mod(Freqs_Pres.DN,1));
    Freq_Scale = sort(unique([Freq_Scale, Freqs_Pres.DN]));
else
    M_aug      = M;
end

        Prob.FL       = Probs_F2L;
        Freqs_Pres.FL = [Rec_Column - Freq_Diff/2; Rec_Column + Freq_Diff/2]; 
    
        Prob.FH       = 1 - Probs_F2L;
        Freqs_Pres.FH = [Rec_Column + Freq_Diff/2; Rec_Column - Freq_Diff/2]; 

switch nev_cond
    case('Low')
        Probs = Prob.L;
        Freqs_Pres = Freqs_Pres.L;
    case('High')
        Probs = Prob.H;
        Freqs_Pres = Freqs_Pres.H;
    case('Equal')
        Probs = Prob.E;
        Freqs_Pres = Freqs_Pres.E;
    case('Diverse Broad')
        Probs = Prob.DB;
        Freqs_Pres = Freqs_Pres.DB;
    case('Diverse Narrow')
        Probs = Prob.DN;
        Freqs_Pres = Freqs_Pres.DN;
    case('Deviant Alone F1')
        Probs = Prob.DA1;
        Freqs_Pres = Freqs_Pres.DA1;
    case('Deviant Alone F2')
        Probs = Prob.DA2;
        Freqs_Pres = Freqs_Pres.DA2;    
    case('Fixed Low')
        Probs = Prob.FL;
        Freqs_Pres = Freqs_Pres.FL;
    case('Fixed High')
        Probs = Prob.FH;
        Freqs_Pres = Freqs_Pres.FH;
end

Freqs_pres_inds = zeros(1,length(Freqs_Pres));

for frin = 1:length(Freqs_Pres)
    Freqs_pres_inds(frin) = find(Freq_Scale == Freqs_Pres(frin));
end

%% Inserting the additional frequencies into h
h_aug = zeros(P,NE,M_aug);
m_count = 1;
for fr = 1:length(Freq_Scale)
    locat = find((1:M) == Freq_Scale(fr));
    if ~size(locat,2)
        h_aug(:,:,fr) = h(:,:,m_count) + (Freq_Scale(fr) - Freq_Scale(fr-1))*(h(:,:,m_count+1)-h(:,:,m_count));
    else
        h_aug(:,:,fr) = h(:,:,m_count);
        m_count = m_count + 1;
    end
end

h = h_aug;

if plot_TC && all_plots
    figure(5)
        h_neur = 67;
        h_col  = 10;
        aug_TC = h(h_neur,h_col,:);
        plot(Freq_Scale,reshape(aug_TC,1,M_aug),'-b.')
        title(['Typical Tuning Curve (Neuron ' num2str(h_neur) ' of Column ' num2str(h_col) ')'])
        xlabel('Frequency')
        ylabel('Input as Fraction of the Input at the Best Frequency')
end

%% Oddball Sequence Generation

switch nev_cond
    case('Low')
        Oddball = [Freqs_Pres(2)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(1)*ones(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('High')
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(2)*ones(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Equal')
        Oddball = [Freqs_Pres(2)*ones(1,ceil(Probs*n_stim)) Freqs_Pres(1)*ones(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Diverse Broad')
        freq_reps = n_stim/length(Freqs_Pres);
        Oddball   = repmat(Freqs_Pres,1,freq_reps);
        Oddball   = Oddball(randperm(n_stim));
    case('Diverse Narrow')
        freq_reps = n_stim/length(Freqs_Pres);
        Oddball   = repmat(Freqs_Pres,1,freq_reps);
        Oddball   = Oddball(randperm(n_stim));
    case('Deviant Alone F1')
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Deviant Alone F2')
        Oddball = [Freqs_Pres(1)*ones(1,ceil(Probs*n_stim)) zeros(1,ceil((1-Probs)*n_stim))];
        Oddball = Oddball(randperm(n_stim));
    case('Fixed Low')
            fix_int = 1/Probs;
            Oddball = Freqs_Pres(1)*ones(1,n_stim);
            Oddball(ceil(rand(1)*fix_int):fix_int:end) = Freqs_Pres(2);
    case('Fixed High')
            fix_int = 1/Probs;
            Oddball = Freqs_Pres(2)*ones(1,n_stim);
            Oddball(ceil(rand(1)*fix_int):fix_int:end) = Freqs_Pres(1);
        
end

Spec_Temp   = zeros(1,num_steps,length(Freqs_pres_inds)); % Spectro-temporal structure of the oddball stimulus sequence
ramp_step   = dt/ramp_dur;
Single_Stim = [ramp_step:ramp_step:1, ones(1,floor((duration - 2*ramp_dur)/dt)), wrev(ramp_step:ramp_step:1)]; % Waveform of each single stimulus
Stim_Onsets = zeros(1,n_stim);
Time_Ind    = zeros(n_stim,2);

for ns = 1:n_stim    
    Time_Ind(ns,:) = floor(t_eq/dt) + floor((ISI+duration)/dt)*(ns - 1) + [0,(length(Single_Stim) - 1)]; % Indices in which there is stimulus
    if Oddball(ns)
    Spec_Temp(1,Time_Ind(ns,1):Time_Ind(ns,2),Freqs_Pres == Oddball(ns)) = Single_Stim;
    end
    Stim_Onsets(ns) = Time_Ind(ns,1);
end

h_sq = h(:,:,Freqs_pres_inds);

%% Returning to equilibrium conditions:
E = E_eq;
I = I_eq;
x = x_eq; 
y = y_eq;
z = z_eq*ones(P,NE,length(Freqs_pres_inds));

%z_I  = z;

s_z  = zeros(size(z));
%s_z_I  = zeros(size(z_I));

Gain_E = Gain_E_eq;
Gain_I = Gain_I_eq;

%E_mean_act = E_mean;
% The following lines are for tracking activity of all neurons (made into comments for memory limitations)
%E_act_overall = zeros(P,NE,num_steps); 
%E_act_overall(zeros(P,NE,M_aug):,:,1:floor(t_eq/dt)) = E_act;

E_mean = [E_mean zeros(P,num_steps - num_steps_eq)];
%I_mean = [I_mean zeros(P,num_steps - num_steps_eq)];

%% Dynamic Loop
if stim
i = floor(t_eq/dt) + 1;
while i < num_steps
    
    % Pre-calculation for the inter-column exc2exc gain calculations:
    EUx    = diag(E*(U.*x)'); % Intra-column input from excitatory synapses, taking synaptic depression into account
    EUx_1 = [EUx(2:P) ; EUx(1)] + [EUx(P); EUx(1:P-1)]; % Excitatory input from neighboring column
    EUx_2 = [EUx(3:P) ; EUx(1:2)] + [EUx(P-1:P); EUx(1:P-2)]; % Excitatory input from the neighboring column once removed

    % Pre-calculation for the intra-column inhibitory gain: 
    IUy = diag(I*(U.*y)'); % Intra-column input from inhibitory synapses, taking synaptic depression into account
    
    % Pre-calculation for the inter-column exc2inh gain:
    E_sum   = sum(E,2); % Total intra-column input in each column
    E_sum_1 = [E_sum(2:P); E_sum(1)] + [E_sum(P); E_sum(1:P-1)]; % Total input to each column from its nearest neighbors
    E_sum_2 = [E_sum(3:P); E_sum(1:2)] + [E_sum(P-1:P); E_sum(1:P-2)]; % Same, for 2nd-nearest neighbors.
     
    % The gain each cell receives:
    Gain_E = bsxfun(@plus,Inp_E,(Jee.*EUx + Jei.*IUy + Jee_1.*EUx_1 + Jee_2.*EUx_2));
    Gain_I = bsxfun(@plus,Inp_I,(Jie.*E_sum + Jii.*sum(I,2) + Jie_1.*E_sum_1 + Jie_2.*E_sum_2));
    %Gain_E = (Jee.*EUx + Jei.*IUy + Jee_1.*EUx_1 + Jee_2.*EUx_2)*ones(1,NE) + Inp_E;
    %Gain_I = (Jie.*E_sum + Jii.*sum(I,2) + Jie_1.*E_sum_1 + Jie_2.*E_sum_2)*ones(1,NI) + Inp_I;
    
    % Calculating the sensory input to all neurons at the current time-step:
    s_E = U_s*A*repmat(Spec_Temp(1,i,:),P,NE).*z.*h_sq;
    
    %s_I = U_s*A*repmat(Spec_Temp(1,i,:),P,NI).*z.*h_sq_I;
    
    % Adding the sensory input to Gain_E:
    Gain_E = Gain_E + stim*sum(s_E,3);
    %Gain_I = Gain_I + inp2inh*stim*sum(s_I,3);
    
    % Dynamics of thalamo-cortical synaptic depression:
    z = z + dt*((1 - z)./tau_rec_s - s_E);
    %z_I = z_I + dt*((1 - z_I)./tau_rec_s - s_I);
    
    % Implementing the non-linearity:
    Gain_E(Gain_E < 0)   = 0;
    Gain_E(Gain_E > 300) = 300;
    Gain_I(Gain_I < 0)   = 0;
    Gain_I(Gain_I > 300) = 300;
    
    % The variables' dynamics:
    E = E + (dt/tau_E)*(-E + Gain_E.*(1 - tau_ref_E*E));
    x = x + dt*((1 - x)./tau_rec - U.*E.*x);         
    I = I + (dt/tau_I)*(-I + Gain_I.*(1 - tau_ref_I*I));
    y = y + dt*((1 - y)./tau_rec - U.*I.*y); 
        
    % Tracking activity (mean and perhaps overall and individual as well):
    E_mean(:,i) = mean(E,2);
    %I_mean(:,i) = mean(I,2);
    
    i = i + 1;
end
end

%% Saving
if save_results    
    
    filename2 = [foldnam '/TYLT_' term '_net' num2str(net_num) '_' nev_cond_code{nevco}];
        
    %if length(Probs_Arr) > 1
    %   filename2 = [filename2 '_P' num2str(Probs_F2L)];
    %end
            
    %if length(ISI_Arr) > 1
    %    filename2 = [filename2 '_ISI' num2str(ISI*1000)];
    %end
    
    %if length(A_Arr) > 1
    %    filename2 = [filename2 '_A' num2str(A)];
    %end

    %filename2 = [filename2 '_tr' num2str(tr) '.mat'];
    filename2 = [filename2 '.mat'];
    save(filename2) % Saving the whole workspace, since certain conditions can lead to some of the variables not being defined.
    
end
curr_time2 = clock;
%disp([num2str(curr_time2(4)) ':' num2str(curr_time2(5)) ' ' nev_cond ' condition, trial ' num2str(tr) ' done']) 
disp([num2str(curr_time2(4)) ':' num2str(curr_time2(5)) ' ' nev_cond ' condition, ISI = ' num2str(ISI) 's, A = ' num2str(A) 'Hz done']) 
toc
E_mean    = [];
I_mean    = [];
Spec_Temp = [];