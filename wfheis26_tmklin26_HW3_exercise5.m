%{
Will Heisler and Tait Kline
HW3 Exercise 5 
11/12/2025
%}

%{ 
       REFERENCE DICTIONARY
C             :   capacitance (uF/cm^2)
I             :   applied current (nA)
I_K           :   potassium channel current (nA)
I_L           :   leakage current (nA)
I_L_inital    :   initial leakage current (nA)
I_Na          :   sodium channel current (nA)
I_P           :   Na-K-ATPase pump current (nA)
V             :   action potential (mV)
V_K           :   displacement from K+ eq potential (mV)
conc_K_in     :   K+ concentration inside (mM/L)
conc_K_out    :   K+ concentration outside (mM/L)
V_Na          :   displacement from Na+ eq potential (mV)
conc_Na_in    :   Na+ concentration inside (mM/L)
conc_Na_out   :   Na+ concentration outside (mM/L)
V_L           :   displacement from leakage eq potential (mV)
g_K           :   maximum K conductance (mS/cm^2)
g_Na          :   maximum Na conductance (mS/cm^2)
g_L           :   maximum leakage conductance (mS/cm^2)
n             :   K activation gating (prob. of K gate being open)
m             :   Na activation gating (prob. of Na gate being open)
h             :   Na inactivation gating (prob. of Na gate being closed)
a_n           :   opening rate constant (ms^-1)
a_m           :   opening rate constant (ms^-1)
a_h           :   opening rate constant (ms^-1)
beta_n           :   closing rate constant (ms^-1)
beta_m           :   closing rate constant (ms^-1)
beta_h           :   closing rate constant (ms^-1)
%}

g_K = 36;
g_Na = 120;
g_L = 0.3;
n = 0.317;
m = 0.05;
h = 0.6; 

V = -65;
V_K = -77;
V_Na = 50;
V_L = -54.4;

C = 0.1;
I = 0;
I_K = g_K*n^4*(V - V_K);
I_Na = g_Na*m^3*h*(V - V_Na);
I_L = g_L*(V - V_L);
I_L_inital = g_L*(V - V_L);

%%% Ion Concentrations %%%
conc_K_in = 150; % mM/L
conc_K_out = 5.5; % mM/L

conc_Na_in = 15; % mM/L
conc_Na_out = 160; % mM/L

% Siumlation parameters %
sim_length = 3; % mS
time_step = 0.001;
num_iter = sim_length / time_step;
time_points = zeros(1, num_iter+1);

% Output data storage %
n_predictions = zeros(1, num_iter+1);
m_predictions = zeros(1, num_iter+1);
h_predictions = zeros(1, num_iter+1);
V_predictions = zeros(1, num_iter+1);

% inital values %
n_predictions(1) = n;
m_predictions(1) = m;
h_predictions(1) = h;
V_predictions(1) = V;

% initialize time %
t = 0;

%%% initialize channel logic %%%
% Activation/deactivation thresholds %
Na_activation = -55;
Na_deactivation = 50; 
K_activation = 50;
L_activation = -54.4; % Leakage channel voltage gating
% Channel open/closed booleans %
Na_open = false;
K_open  = false;
L_open = false; 


for iter = 1:num_iter

    % Applied current time logic %
    if t > 0.5 && t < 1.0 % 
        I = 15;
    else 
        I = 0;
    end

    % Opening rate constants %
    a_n = 0.01 * ((V +55)/(1-exp(-(V+55)/10)));
    a_m = 0.1 * ((V + 40)/(1-exp(-(V+40)/10)));
    a_h = 0.07 * exp(-(V+65)/20);
   
    % Closing rate constants % 
    beta_n = 0.125 * exp(-(V+65)/80);
    beta_m = 4 * exp(-(V+65)/18);
    beta_h = 1/(exp(-(V+35)/10)+1);

    % Gate opening probabilities % 
    dn_dt = @(t, a_n, n, beta_n) a_n*(1-n)-(beta_n*n);
    dm_dt = @(t, a_m, m, beta_m) a_m*(1-m)-(beta_m*m);
    dh_dt = @(t, a_h, h, beta_h) a_h*(1-h)-(beta_h*h);

    %%% voltage-gated logic %%%

    % Close K when V reaches minimum value and begins to increase
     if K_open 
        if V_predictions(iter - 1) < V
            K_open = false;
        end
     end
    % Open K when V reaches 49.3
    if V >= K_activation
        K_open = true;
    end
    % Open Na when V exceedes -55
    if V > Na_activation && ~K_open % only if K isn't already open
        Na_open = true;
    end
    % Close Na when V reaches 49.3
    if V >= Na_deactivation
        Na_open = false;
    end
    % Open leakage channel when V <= -54.4
    if V <= L_activation
        L_open = true;
    else
        L_open = false;
    end

    % Na and K channel current calculations %
    if Na_open
        I_Na = g_Na*(m^3)* h*(V - V_Na);
    else
        I_Na = 0;
    end
    conc_Na_out = conc_Na_out + (I_Na*time_step); 
    conc_Na_in = conc_Na_in - (I_Na*time_step);
    
    if K_open
        I_K = g_K * (n^4) * (V - V_K);
    else
        I_K = 0;
    end
    conc_K_out = conc_K_out + (I_K*time_step);
    conc_K_in = conc_K_in - (I_K*time_step);

    %%% Leakage channel current calculation %%%
    if L_open
        I_L = g_L*(V - V_L);  
    else
        I_L = 0;
    end
    conc_K_out = conc_K_out + (I_L*time_step);
    conc_K_in = conc_K_in - (I_L*time_step); 

    %%% ATP-ase pump current calculation %%%
    %Ion concentration-based pump activation %
    if conc_K_out > 0 && conc_Na_in > 0
        I_P = -1 * I_L_inital;
    else
        I_P = 0;
    end
    conc_Na_out = conc_Na_out + (I_P*3*time_step);
    conc_Na_in = conc_Na_in - (I_P*3*time_step);
    conc_K_out = conc_K_out - (2*I_P*time_step);
    conc_K_in = conc_K_in + (2*I_P*time_step);
    
    % Constituative pump activation %
    %{
    I_P = -1 * I_L_inital;
    conc_Na_out = conc_Na_out + (I_P*3*time_step);
    conc_Na_in = conc_Na_in - (I_P*3*time_step);
    conc_K_out = conc_K_out - (2*I_P*time_step);
    conc_K_in = conc_K_in + (2*I_P*time_step);
    %}
    
    %%% Membrane potential diff. eq. %%%
    dV_dt = @(t, I, I_K, I_Na, I_L, I_P, C) (I - I_K - I_Na - I_L - I_P)/C;

    %%% RK4 Calculations %%%

    % RK4 Step 1
    dn1 = dn_dt(t, a_n, n, beta_n) * time_step;
    dm1 = dm_dt(t, a_m, m, beta_m) * time_step;
    dh1 = dh_dt(t, a_h, h, beta_h) * time_step;
    dV1 = dV_dt(t, I, I_K, I_Na, I_L, I_P, C) * time_step;
    
    % RK4 Step 2
    dn2 = dn_dt(t + 0.5*time_step, a_n, n + 0.5*dn1, beta_n) * time_step;
    dm2 = dm_dt(t + 0.5*time_step, a_m, m + 0.5*dm1, beta_m) * time_step;
    dh2 = dh_dt(t + 0.5*time_step, a_h, h + 0.5*dh1, beta_h) * time_step;
    dV2 = dV_dt(t + 0.5*time_step, I, I_K, I_Na, I_L, I_P, C) * time_step;

    % RK4 Step 3
    dn3 = dn_dt(t + 0.5*time_step, a_n, n + 0.5*dn2, beta_n) * time_step;
    dm3 = dm_dt(t + 0.5*time_step, a_m, m + 0.5*dm2, beta_m) * time_step;
    dh3 = dh_dt(t + 0.5*time_step, a_h, h + 0.5*dh2, beta_h) * time_step;
    dV3 = dV_dt(t + 0.5*time_step, I, I_K, I_Na, I_L, I_P, C) * time_step;
    
    % RK4 Step 4
    dn4 = dn_dt(t + time_step, a_n, n + dn3, beta_n) * time_step;
    dm4 = dm_dt(t + time_step, a_m, m + dm3, beta_m) * time_step;
    dh4 = dh_dt(t + time_step, a_h, h + dh3, beta_h) * time_step;
    dV4 = dV_dt(t + time_step, I, I_K, I_Na, I_L, I_P, C) * time_step;
    
    % Change in variable 
    delta_n = (1/6) * (dn1 + 2*dn2 + 2*dn3 + dn4);
    delta_m = (1/6) * (dm1 + 2*dm2 + 2*dm3 + dm4);
    delta_h = (1/6) * (dh1 + 2*dh2 + 2*dh3 + dh4);
    delta_V = (1/6) * (dV1 + 2*dV2 + 2*dV3 + dV4);

    % Update variables
    n = n + delta_n;
    m = m + delta_m;
    h = h + delta_h;
    V = V + delta_V;

    % Update time
    t = t + time_step;

    % Store predicted values
    time_points(iter+1) = iter * time_step; 
    n_predictions(iter+1) = n;
    m_predictions(iter+1) = m;
    h_predictions(iter+1) = h;
    V_predictions(iter+1) = V;

    %%% Concentration sampling %%%
    %{
    if abs(t- 0.001) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);  
    end 

    if abs(t - 0.4) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);
    end

    if abs(t - 0.6) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);
    end

    if abs(t - 0.8) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);
    end

    if abs(t - 1.0) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);
    end

    if abs(t - 2.63) < 0.0001
        fprintf('K_in: %.3f, K_out: %.3f, Na_in: %.3f, Na_out: %.3f\n', ...
        conc_K_in, conc_K_out, conc_Na_in, conc_Na_out);
    end
    %}

end
    %%% Figures %%%
    
    % Gate opening probability plot %
    figure;
    plot(time_points, n_predictions,'r.' , 'LineWidth',1.5)
    hold on
    plot(time_points, m_predictions,'b.' , 'LineWidth',1.5)
    hold on 
    plot(time_points, h_predictions,'g.' , 'LineWidth',1.5)
    legend({'n (prob. of K gate being open)', 'm (prob. of Na gate being open)',...
        'h (prob. Na gate being closed)'}, 'Location', 'best');
    title('Voltage-Gated Channel Opening Probabilities')
    xlabel('Time (ms)')

    % Membrane potential plot % 
    figure
    plot(time_points, V_predictions,'r.' , 'LineWidth',1.5)
    title('Neuron Membrane Potential Over Time')
    ylabel('Membrane potential (mV)')
    xlabel('Time (ms)')