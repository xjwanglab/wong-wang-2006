%% Code for simulating a reduced mean-field model for decision-making (see
%% Appendix of Wong & Wang, J. Neurosci. (2006))

%%%% Synaptic time and other constants %%%%%%%%%%%%%%%%%%

Tnmda = 100;    % NMDAr
Tampa = 2;      % AMPAr
gamma = 0.641;

%%%% FI curve parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 270; b = 108; 
d = 0.1540;  % Parameters for excitatory cells

%%%% Vectorizing variables for evaluations at end of (block) loop %%%%%%%%%%

r1_traj = [];  r2_traj = []; 
s1_traj = [];  s2_traj = []; 
cross_r1 = []; cross_r2 = [];

%%%% Parameters to be varied %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh       = 6.4;       % Coherence level
mu0       = 30.0;      % External stimulus strength
thresh    = 15;        % Decision threshold
thresh_s  = (gamma*Tnmda*thresh/1000)/(1+gamma*Tnmda*thresh/1000); % Threshold in s-space
noise_amp = 0.02;      % Noise amplitude into selective populations
N_trials  = 50 ;       % Total number of trials
Tstim     = 1500;      % Stimulus duration (ms)

%%%%% Stimulus input strengths %%%%%

mu1 = mu0*(1+coh/100); % input strength to pop 1 (for coherence)
mu2 = mu0*(1-coh/100); % input strength to pop 2 (for coherence)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Trials number and (block) loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ww = 1:N_trials % Trial loop

    trial_no = ww; 

    %---- Vectorise variables with sliding window -------------
    nu1_wind = [] ; nu2_wind = [] ; 
    s1_wind  = [] ; s2_wind  = [] ; 

    %---- Initial conditions and clearing variables -----------
    s1_in=0.1; s2_in=0.1;
    clear nu1_in nu2_in I_eta1 I_eta2  ;
    nu1_in = 2; nu2_in = 2; 
    I_eta1_in = noise_amp*randn ; I_eta2_in = noise_amp*randn ;

    %---- Time conditions -------------------------------------
    
    dt = 0.5;           % Time step in msec
    T_total = 3000/dt;  % Total number of steps
    time_wind = 50/dt;  % Temporal window size for averaging
    slide_wind = 5/dt;  % Sliding step for window

    %---- Intialise and vectorise variables to be used in loops below ------
    
    s1 = s1_in.*ones(1,T_total); s2 = s2_in.*ones(1,T_total); 
    nu1 = nu1_in.*ones(1,T_total); nu2 = nu2_in.*ones(1,T_total);
    phi1 = nu1_in.*ones(1,T_total); phi2 = nu2_in.*ones(1,T_total); 
    I_eta1 = I_eta1_in.*ones(1,T_total); I_eta2 = I_eta2_in.*ones(1,T_total); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    for t = 1:T_total

        %---- Constant effective external current input (with inhibition taken into account)
        I0E1 = 0.3255; I0E2 = 0.3255;
        
        %---- External stimulus---------------------------------------------------------
        JAext = 0.00052; % Synaptic coupling constant to external inputs
        I_stim_1 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu1); % To population 1
        I_stim_2 = (500/dt<t & t<(Tstim+500)/dt)*(JAext*mu2); % To population 2

        %---- Recurrent synaptic coupling constants-------------------------------------
        JN11 = 0.2609; JN22 = 0.2609; 
        JN12 = 0.0497; JN21 = 0.0497; 
        
        %---- Resonse function of competiting excitatory population 1 ------
	    Isyn1(t) = JN11.*s1(t) - JN12.*s2(t) + I0E1 + I_stim_1 + I_eta1(t) ;
	    phi1(t)  = (a.*Isyn1(t)-b)./(1-exp(-d.*(a.*Isyn1(t)-b)));

	    %---- Response function of competiting excitatory population 2 -----
	    Isyn2(t) = JN22.*s2(t) - JN21.*s1(t) + I0E2 + I_stim_2 + I_eta2(t) ;
	    phi2(t)  = (a.*Isyn2(t)-b)./(1-exp(-d.*(a.*Isyn2(t)-b)));

	    %---- Dynamical equations -------------------------------------------

	    % Mean NMDA-receptor dynamics
	    s1(t+1) = s1(t) + dt*(-(s1(t)/Tnmda) + (1-s1(t))*gamma*nu1(t)/1000);
	    s2(t+1) = s2(t) + dt*(-(s2(t)/Tnmda) + (1-s2(t))*gamma*nu2(t)/1000);

        % Noise through synaptic currents of pop1 and 2
        I_eta1(t+1) = I_eta1(t) + (dt/Tampa)*(-I_eta1(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        I_eta2(t+1) = I_eta2(t) + (dt/Tampa)*(-I_eta2(t)) + sqrt(dt/Tampa)*noise_amp*randn ;

        % To ensure firing rates are always positive (noise may cause negative)
        if phi1(t) < 0 
            nu1(t+1) = 0;
            phi1(t) = 0;
        else
            nu1(t+1) = phi1(t);
        end;
        if phi2(t) < 0
            nu2(t+1) = 0;
            phi2(t) = 0;
        else
            nu2(t+1) = phi2(t);
        end;
        
	    %==============================================================================================
        
    end;  %---- End of time loop --------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %---- Calculating the mean rates and gating variables with sliding window -----

    nu1_wind = [nu1_wind (mean(nu1(1:time_wind)))] ;
    nu2_wind = [nu2_wind (mean(nu2(1:time_wind)))] ;
    s1_wind  = [s1_wind (mean(s1(1:time_wind)))] ;
    s2_wind  = [s2_wind (mean(s2(1:time_wind)))] ;

    for t = 1:(T_total-time_wind)/slide_wind

        nu1_wind = [nu1_wind (mean(nu1(slide_wind*t:slide_wind*t+time_wind)))] ;
        nu2_wind = [nu2_wind (mean(nu2(slide_wind*t:slide_wind*t+time_wind)))] ;
        s1_wind  = [s1_wind (mean(s1(slide_wind*t:slide_wind*t+time_wind)))] ;
        s2_wind  = [s2_wind (mean(s2(slide_wind*t:slide_wind*t+time_wind)))] ;
        
    end;

    r1_traj=[r1_traj; nu1_wind]; r2_traj=[r2_traj; nu2_wind];
    s1_traj=[s1_traj; s1_wind]; s2_traj=[s2_traj; s2_wind];
    clear nu1 nu2 s1 s2; 
    clear nu1_wind nu2_wind s1_wind s2_wind ;

end; %---- End trial loop ---------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Plots of first few trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_traj = 10; % Total number of trajectories

subplot(2,2,1:2) % Plot timecourse of activity (population firing rates)
for ww = 1:N_traj
    plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r1_traj(ww,1:end-10),'b');hold on;
    plot([dt*time_wind:dt*slide_wind:dt*(T_total-time_wind)],r2_traj(ww,1:end-10),'r--');hold on;
end;
%grid on; 
%axis tight;
hold on; plot([1 dt*(T_total-time_wind)],[thresh thresh],'k--');
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
title('Timecourse of firing rates r_1 (blue) and r_2 (red)'); 

subplot(2,2,3) % Plot trajectories in phase-space
for ww=1:N_traj
    if s1_traj(ww,end)>thresh_s
       plot(s2_traj(ww,:),s1_traj(ww,:),'b');hold on;
    else
       plot(s2_traj(ww,:),s1_traj(ww,:),'r');hold on;
    end;
end;
hold on; plot([0 1],[thresh_s thresh_s],'k--');
hold on; plot([thresh_s thresh_s],[0 1],'k--');
xlabel('S_2'); ylabel('S_1'); 

subplot(2,2,4) % Plot trajectories in phase-space
for ww=1:N_traj
    if r1_traj(ww,end)>thresh
       plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'b');hold on;
    else
       plot(r2_traj(ww,1:end),r1_traj(ww,1:end),'r');hold on;
    end;
end;
hold on; plot([0 45],[thresh thresh],'k--');
hold on; plot([thresh thresh],[0 45],'k--');
xlabel('r_2 (Hz)'); ylabel('r_1 (Hz)'); 
axis tight;
