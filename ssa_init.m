%%  Application of Loebel's model for SSA
% The following code has been adapted from that of Alex Loebel
% Scalars are denoted by letters as in the article, or else by lower-case words.
% Vectros and matrices are denoted by letters as in the article, or else by capitalized words.

% This file generates and brings to equilibrium a network that can be used
% to run Nevo's protocols

%close all
%clear all

all_plots   = 0; % Controls which plots to display
plot_mixing = 1; % Plots the best frequencies for all neurons
plot_TC     = 1; % Plots the shape of a typical tuning curve

find_FRA = 0;
id_crit = 0;

Sel_Columns = [8 9 10]; % Selects columns whose mean activity will be plotted
n_sing      = 1; % Number of single neurons that will be tracked in each selected column
inc_non_act = 1; % Includes non-active neurons in the single neurons tracked

graph_col = 'bgrcmyk'; % Colors for the graphs of the different cases
n_sel     = length(Sel_Columns);

Conds{1} = 'Low';
Conds{2} = 'High';
Conds{3} = 'Equal';
Conds{4} = 'Diverse Broad';
Conds{5} = 'Diverse Narrow';
Conds{6} = 'Deviant Alone F1';
Conds{7} = 'Deviant Alone F2';
Conds{8} = 'Fixed Low';
Conds{9} = 'Fixed High';

nev_cond_code{1} = 'L';
nev_cond_code{2} = 'H';
nev_cond_code{3} = 'E';
nev_cond_code{4} = 'DB';
nev_cond_code{5} = 'DN';
nev_cond_code{6} = 'DA1';
nev_cond_code{7} = 'DA2';
nev_cond_code{8} = 'FL';
nev_cond_code{9} = 'FH';

%% Parameters of the Model:
% Number of columns:
P = 21;

% Number of frequencies to which the network is sensitive
freq_stretch = 1; % Choose an integer >=1; this multiplies P to give M.
M = freq_stretch*P; % Thus, column best frequencies will be certain presentable frequencies.

ring_net = 0;

% Numbers of cells in each column (excitatory and inhibitory):
NE = 100;
NI = 100; 

% Controlling connection strengths:
factor   = 1; % multiplies all intra-column connections (controls spontaneous population spikes)
factor_1 = 1; % multiplies nearest-neighbor inter-column connections (aids spread of population spikes)
factor_2 = 1; % multiplies 2nd-nearest-neighbor inter-column connections (aids spread of population spikes)

% Background input drawn from a uniform distribution (neurons are indexed by input strength):
bg_low_E  = -9.9; % lowest background input (in Hz) to excitatory population
bg_high_E = 9.9; % highest background input (in Hz) to excitatory population
bg_low_I  = bg_low_E; % lowest background input (in Hz) to inhibitory population
bg_high_I = bg_high_E; % highest background input (in Hz) to inhibitory population
% *) In Loebel & Tsodyks 2002, the above values were reached by requiring a
% spontaneous activity of a few Hz. 

% Time constants:
tau_E     = 0.001; % excitatory neurons' time constant (in seconds)
tau_I     = 0.001; % inhibitory neurons' time constant (in seconds) 
tau_ref_E = 0.003; % tau refractory of excitatory neurons (in seconds)
tau_ref_I = 0.003; % tau refractory of inhibitory neurons (in seconds)
tau_rec   = 0.800; % recovery time constant of intracortical synapses (in seconds)
tau_rec_s = 0.300; % recovery time constant of sensory input synapses (in seconds)
% *) This value was chosen since Eli says that LFP recovers fully after about 1 s (tau_rec_s ~ 0.3 s); there should be
% articles about measurements in slices showing this recovery

U   = 0.5; % Portion of available fraction of resources that is utilized in response to an action potential
U_s = 0.7; % Same as U, only for the thalamo-cotical synapses that convey the sensory input

% Connection strengths:
Jee = 6*factor/NE;    % exc2exc
Jei = -4*factor/NI;   % inh2exc
Jie = 0.5*factor/NE;  % exc2inh
Jii = -0.5*factor/NI; % inh2inh

Jee_1 = 0.045*factor_1/NE;  % exc2exc, between neighboring columns
Jie_1 = 0.0035*factor_1/NE; % exc2inh, between neighboring columns
Jee_2 = 0.015*factor_2/NE;  % exc2exc, from one column to its 2nd neighbor
Jie_2 = 0.0015*factor_2/NE; % exc2inh, from one column to its 2nd neighbor

% Simulation conditions:
t_eq      = 5; % The time given to reach equilibrium (in seconds). It is important to allow enough time, or else the response to the first stimulus is distorted.
dt        = 0.0001; % Time-step (in seconds)
post_stim = 0.100; % This is defined here so the simulations keeps runnng for 2*post_stim after the last stimulus offset

%% Stimulus Parameters
A = 5; % Peak magnitude of the input (the magnitude felt by neurons most sensitive to that input)
A_max = 50; % This is used for tuning curves at high amplitudes (comes into effect if A_max > alpha)

ISI  = 0.3; % Inter-stimulus interval (in seconds)
duration = 0.050; % Total duration of each stimulus (in seconds)
ramp_dur = 0.005; % durations of the ramps at the beginning and end of each stimulus (in seconds)

n_stim = 100; % Total no. of stimuli (Best take a product of 10)
%t_prot = n_stim*(duration + ISI); % Total time of the oddball protocol

%% Localization of Sensory Input Effect:
lambda_c = 5; %0.25; % Base value of the lambda_s localization parameter
alpha    = 101; %2; % Magnitude threshold for increase in localization parameter
delta    = 5; % Slope of the increase in lambda_s for magnitudes greater than alpha (determines the degree of localization of an input at high sound levels)

% Option for asymmetric tuning curves:
delta_right = delta; 
delta_left  = delta;

%% Mixing cells with different frequency preferences in each column:
mix_tcs = 1; % 1 or 0, mixing or segregating the input to individual columns

part_2f = 16; % 1 over the fraction of cells in each column that are sensitive to the frequency of the column 2 columns forward is about 1/12
part_2b = 16; % Same, but for the column 2 columns back
part_1f = 8; % 1 over the fraction of cells in each column that are sensitive to the frequency of the next column is about 1/6
part_1b = 8; % Same, but for 1 column back

%%
act_thresh     = 0; % Defines a neuron as active or non-active. In Loebel et al. 2007, this was chosen as 0.
%% 
cancel_syns = 0;
pr_syn_canc = 0.5; % Probability of synapse cancelation
overlap     = 0; % Whether or not there may be neurons responsive to both F2 and F1
zipper      = 1; % This means that canceled synapses in the specified columns corresponding to stimulus frequencies follow a zipper pattern
first_canc  = 1; % In case of zipper pattern, the first synapse canceled in F1
        
%% Sensory Input
lambda_c = 2.5/(log(2)^0.5);

lambda_s_right = lambda_c + (A_max > alpha)*(A_max - alpha)/delta_right; % Determines localiztion of input effect over the cortical sheet
lambda_s_left  = lambda_c + (A_max > alpha)*(A_max - alpha)/delta_left;

if ring_net
    Distances = toeplitz([0:floor(M/2) wrev(1:(floor(M-1)/2))]); % Distance of each frequency from all the others.
else
    Distances = toeplitz(0:(M-1));
end

Right_Dec_Args = Distances./lambda_s_right;
Left_Dec_Args  = Distances./lambda_s_left;

curve_type = 'linear'; % Type of the tuning curve; see cases for options and details.

switch(curve_type)
    case 'exponential' % this one is the original type, used in Loebel et al. 2007. 
        h_outline = (triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)); % Magnitudes of input received by each column for each presented frequency
    
    case 'biased_exponential' % this is similar, but would be sharper and have inhibition from all farther off frequencies
        curve_bias = 1/3;
        h_outline  = (1 + curve_bias)*(triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)) - curve_bias; % Ensures a magnitude of A at the tuning curve's peak
    
    case 'linear' % This type is triangular, following Eli's request that some columns do not receive input at all.
        h_outline = (triu(1-Right_Dec_Args) + tril(1-Left_Dec_Args,-1));
        h_outline(h_outline < 0) = 0;
        
    case 'mexican hat' % My own addition; has support in STRFs such as in deCharms et al. 1998; Loebel et al. 2007 state that lateral supression arises in the model due to synaptic depression 
        mex_dec = 1;
        mex_fre = 1;
        h_outline = (triu(cos(mex_fre*Right_Dec_Args).*exp(-mex_dec*Right_Dec_Args)) + tril(cos(mex_fre*Left_Dec_Args).*exp(-mex_dec*Left_Dec_Args),-1));

    case 'Gaussian' 
        h_outline = (triu(exp(-Right_Dec_Args.^2)) + tril(exp(-Left_Dec_Args.^2),-1)); % Magnitudes of input received by each column for each presented frequency
    
end

h_outline = h_outline(1:freq_stretch:end,:); % This gives the columns best frequencies that are spread evenly across the range.
h_outline(P+1:end,:) = []; % Cropping the rows of h_outline to the no. of columns, just in case.

h_outline = reshape(h_outline,[P 1 M]); % *) Aligning the tuning curves to the neurons.
h = repmat(h_outline,[1 NE 1]);

inp2inh = 0; % This states whether there is sensory input to inhibitory neurons

%% Plotting the Tuning Curve 
if plot_TC %&& all_plots
    tc_column = Sel_Columns(3); %floor(P/2);
    
    figure(1)
        hold on
        plot(reshape(h_outline(tc_column,1,:),1,M),'-b.')

        title(['Typical Tuning Curve (Column ' num2str(tc_column) ')'])
        xlabel('Frequency')
        ylabel('Input as Fraction of the Input at the Best Frequency')
        xlim([0 22])
end

%% Mixing the Responses of Neurons in Each Column
if mix_tcs
   h_unmixed = h; % Original h is saved for FRA finding in unmixed case.
   shif = freq_stretch; % The basic shift in BF.
   placer1 = floor(NE/part_1f);
   placer2 = placer1 + floor(NE/part_1b);
   placer3 = placer2 + floor(NE/part_2f);
   placer4 = placer3 + floor(NE/part_2b);
   
   for Q = 1:P % Doing an independent scramble for each column
        Mixed   = randperm(NE);
        
        h(Q,Mixed(1:placer1),:)              = h(Q,Mixed(1:placer1),[((M-shif+1):M) 1:(M-shif)]); % forward means to the right, i.e. increase in 3rd dimension
        h(Q,Mixed(1:placer1),((M-shif+1):M)) = ring_net*h(Q,Mixed(1:placer1),((M-shif+1):M));
        
        h(Q,Mixed((placer1+1):placer2),:)      = h(Q,Mixed((placer1+1):placer2),[(shif+1):M (1:shif)]);
        h(Q,Mixed((placer1+1):placer2),1:shif) = ri
  