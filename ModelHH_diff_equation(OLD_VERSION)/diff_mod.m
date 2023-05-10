%% Conversion to SI
global cm mV ms mS uF uA
    cm = 1e-2;
    mV = 1e-3;
    ms = 1e-3;
    mS = 1e-3;
    uF = 1e-6;
    uA = 1e-6;

%% General value for the system
global gNa eNa gK eK gL eL Cm startTime Ipulse Tpulse I0 Vm_rest Cm_t_out count z start_time_Delta_T rise_time Delta_Temp tau
    gNa = 120 * mS/(cm*cm);  
    eNa = -115 * mV; %Change the sign of the potentiel because of the convention of the model
    gK = 36 * mS/(cm*cm);  
    eK = 12 * mV;
    gL = 0.3 * mS/(cm*cm);  
    eL = -10.6 * mV;
    Cm = 1 * uF/(cm*cm); 
    startTime = 5*ms;
    Ipulse = 75*(uA/(cm*cm));
    Tpulse = 2*ms;
    I0 =0*uA/(cm*cm); 
    Vm_rest= -65*mV;
    Cm_t_out = [];
    count = 0;
    start_time_Delta_T = 10*ms;
    rise_time = 1*ms;
    Delta_Temp = 25;
    tau = 0.1* 1e2;



SimT = 80*ms;
odeTime = 0:0.01*ms:SimT;       % Simulation time steps

VmCI = 0*mV;
nCI=0.3177;                     % These values are obtained with the previous graph
mCI=0.05296;                    % for Vm=0 because it's the initial condition stated above. 
hCI=0.596;
CmCI = Cm;
Cm_t = [];

odeCI = [VmCI nCI mCI hCI];


%Parameters for the Temperature pulse



Temp_init = 20;

z = test_model_T(Delta_Temp, rise_time*1000, odeTime, start_time_Delta_T*1000, Temp_init);
disp(length(z));
figure
plot(z);

 

[t,y]=ode15s(@dy, odeTime, odeCI);
%[t,y] = ode15s(@(t,y) dy(t,y, Cm_t), odeTime, odeCI);
Vm = -y(:,1);
n4 = y(:,2).^4;
m3 = y(:,3).^3;
h = y(:,4);

figure
plot(Cm_t_out);


stimulationCurrent = generateStimCurrent(t*1000,startTime*1000, Tpulse*1000, Ipulse, I0); %generates current, t converted in second because defined in ms
% chartitle=sprintf('Time independent current: I0 = %d µA/cm²',I0/(uA/(cm*cm)));


Vm = (Vm*1000) - 65;
t = t*1000;
figure
plot(t, Vm);
ylabel("Potential membrane (mV)");
xlabel("Time (ms)");
chartitle=sprintf('Pulse current: Ipulse = %d µA/cm², Tpulse = %d ms at starting time = %d ms', Ipulse/(uA/(cm*cm)), Tpulse/ms, startTime/ms);
title(chartitle);
% hold on
% plot(t, stimulationCurrent);




% Function to implemente the differential equation
function dy = dy(t, y)
    global Cm startTime Tpulse Ipulse I0 Cm_t_out Vm_rest
    dy = zeros(4,1);

    disp(t);

    Iinj = generateStimCurrent(t*1000, startTime*1000, Tpulse*1000, Ipulse, I0);
    %Cinj = Cm_inj((t*1000)+1);

    Vm = y(1);
    n =  y(2);
    m =  y(3);
    h =  y(4);

    
    %Cm_t_out= [Cm_t_out, Cinj];

    

    dy(1) = (1/Cm)*( -Iinj - INa(Vm,m,h) - IK(Vm, n) - IL(Vm));                  %dV/dt
    %dy(1) = (1/Cm)*( -Iinj - INa(Vm,m,h) - IK(Vm, n) - IL(Vm)) - dCmdt((t*1000)+1); 
    dy(2) = (alphaN(Vm)*(1-n) - betaN(Vm)*n)*1000;                                                        %dn

    dy(3) = (alphaM(Vm)*(1-m) - betaM(Vm)*m)*1000;                                                        %dm 

    dy(4) = (alphaH(Vm)*(1-h) - betaH(Vm)*h)*1000;                                                        %dh 
end




%Function for the equation of the capacitance
function Cminj = Cm_inj(t)
    global uF cm z
    Cminj = 0.824 * uF/(cm*cm) + ((2.2 * uF/(cm*cm))/(31-z(t)));
end

function dCmdt = dCmdt(t)
    global uF z start_time_Delta_T Delta_Temp rise_time tau
    start_time_Delta_T = start_time_Delta_T*1000;
    rise_time = rise_time*1000;

    if t < start_time_Delta_T
        dCmdt = 0;
    end

    if t > start_time_Delta_T && t < rise_time
        pente = Delta_Temp/rise_time;
        dCmdt = (2.2 * uF*pente) / ((31-z(t))^2);
    end

    if t > rise_time
        dCmdt = (2.2 * uF* Delta_Temp * (-1/tau) * exp(-((t-rise_time-start_time_Delta_T)/tau)))/(31-z(t))^2;
    end
end


%Function for the different ions current
function INa = INa(V, m, h)
    global gNa eNa Vm_rest
    %V =  -V + (Vm_rest);
    INa = gNa * m^3 * h * (V - eNa);
end

function IK = IK(V, n)
    global gK eK Vm_rest
    %V = - V + (Vm_rest);
    IK = gK * n^4 * (V - eK); 
end

function IL = IL(V)
    global gL eL Vm_rest
    %V = - V + (Vm_rest);
    IL = gL * (V - eL);
end



% alpha and beta functions for the gating variables 
function aM = alphaM(V)
    global Vm_rest
    %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    aM = (0.1*(25+V))/(exp((25+V)/10)-1);
end

function bM = betaM(V)
    global Vm_rest
      %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    bM = 4*exp(V/18);
end

function aH = alphaH(V)
    global Vm_rest
     %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    aH = 0.07*exp(V/20);
end

function bH = betaH(V)
    global Vm_rest
    %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    bH = 1/(exp((30+V)/10)+1);
end

function aN = alphaN(V)
    global Vm_rest
      %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    aN = 0.01*(10+V)/(exp((10+V)/10)-1); 
end

function bN = betaN(V)
    global Vm_rest
     %V = ((Vm_rest - V)*1000); 
    V = V * 1000;
    bN = 0.125*exp(V/80);
end
