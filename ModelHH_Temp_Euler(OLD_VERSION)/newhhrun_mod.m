function [V,m,h,n,t] = newhhrun_mod(I, tspan, v, mi, hi, ni, Plot)

    %% Conversion to SI
    global cm mV ms mS uF uA
    cm = 1e-2;
    mV = 1e-3;
    ms = 1e-3;
    mS = 1e-3;
    uF = 1e-6;
    uA = 1e-6;

  dt = 0.001;               % time step for forward euler method
  loop  = ceil(tspan/dt);   % no. of iterations of euler
  
  gNa = 120 * mS/(cm*cm);  
  eNa = 115 * mV;
  gK = 36 * mS/(cm*cm);  
  eK = -12 * mV;
  gL = 0.3 * mS/(cm*cm);  
  eL = 10.6 * mV;

  % D = 5* 1e-6*1e-4; % 1 à 15 um
  % d = 0.6*D; % Axon diameter
  % l = 1.5*1e-6*1e-4; % nodal width (in cm)


  % Initializing variable vectors
  t = (1:loop)*dt;
  V = zeros(loop,1);
  m = zeros(loop,1);
  h = zeros(loop,1);
  n = zeros(loop,1);
  Cm = zeros(loop,1);
  INa = zeros(loop,1);
  IK = zeros(loop,1);
  IL = zeros(loop,1);
  Itot = zeros(loop,1);
  IC = zeros(loop,1);

   % Set initial values for the variables
  
  V(1) = v * mV;
  m(1) = mi;
  h(1) = hi;
  n(1) = ni;

  Cm(1) = 1 * uF/(cm*cm);
  T_init = 18.5;


  % Vecteur temperature
  t2 = 0:0.001:1000;
  z = test_model_T(20, 100, t2, 100000, T_init); %Delta_T / Temps de montée / vecteur temps / start time / Temp initial


  % Euler method %%%%%%%ATTENTION, cette méthode fonctionne seulement si
  % les Paramètres n,m,h et les courants sont constant! ICI ce n'est pas le
  % cas à cause de la variation de la température
  % for i=1:loop-1     
  %     Cm(i) = Cm(1) * ((1 + 0.0048*(z(i) - T_init))/(1 - 0.002*(z(i) - T_init))); %Equation source Optocapacitance 
  % 
  %     %Cm(i) = 0.824 * uF/(cm*cm) + ((2.2 * uF/(cm*cm))/(31-z(i)));
  % 
  %     INa(i) = gNa*((m(i))^3)*h(i)*(V(i) - eNa); 
  %     IK(i) = gK*((n(i))^4)*(V(i) - eK); 
  %     IL(i) = gL*(V(i) - eL);
  % 
  %     Itot(i) = (I - (INa(i) + IK(i) + IL(i))) * (uA/(cm*cm));
  % 
  %     V(i+1) = V(i) + ((Itot(i)/Cm(i))*dt)*mV;
  %     m(i+1) = m(i) + dt*(alphaM(V(i) )*(1-m(i)) - betaM(V(i) )*m(i));
  %     h(i+1) = h(i) + dt*(alphaH(V(i) )*(1-h(i)) - betaH(V(i) )*h(i));
  %     n(i+1) = n(i) + dt*(alphaN(V(i) )*(1-n(i)) - betaN(V(i) )*n(i));
  % 
  %     IC(i) = Cm(i) * ((V(i+1) - V(i))/dt);
  % end

  

  l = getpeaks(V);
  disp(l)
  blb = 200000;
  %qs = filt_peak(l,blb);

  if Plot == 1
    V = V*1000;
    figure
    plot(t,V);
    findpeaks(V)
    xlabel('Time');
    ylabel('Membrane Potential');
    title('Voltage time series');
    
    
    figure
    plot(z);
    xlabel('Time');
    ylabel('Temperature');
    title('Temperature vs. time');

    figure
    plot(Cm);
    xlabel('Time');
    ylabel('Capacitance');
    title('Capacitance vs. time');

    figure
    plot(V);
    xlabel('Time');
    ylabel('Voltage membrane');
    title('Voltage membrane vs. time');

    figure
    plot(INa);
    xlabel('Time');
    ylabel('Current channels');
    title('Current channels vs. time');
    hold on;
    plot(IK);
    hold on;
    plot(IL);
    legend("INa", "IK", "IL");

    figure
    plot(Itot);
    xlabel('Time');
    ylabel('Current channels total');
    title('Current channels total vs. time');

    figure
    plot(IC);
    xlabel('Time');
    ylabel('Current in capacitance');
    title('Current in capacitance vs. time');
  end
end

% alpha and beta functions for the gating variables 

function aM = alphaM(V)
V = (-65*1e-3) - V; 
aM = (0.1*(25+V))/(exp((25+V)/10)-1);
end

function bM = betaM(V)
V = (-65*1e-3) - V; 
bM = 4*exp(V/18);
end

function aH = alphaH(V)
V = (-65*1e-3) - V; 
aH = 0.07*exp(V/20);
end

function bH = betaH(V)
V = (-65*1e-3) - V; 
bH = 1/(exp((30+V)/10)+1);
end

function aN = alphaN(V)
V = (-65*1e-3) - V; 
aN = 0.01*(10+V)/(exp((10+V)/10)-1); 
end

function bN = betaN(V)
V = (-65*1e-3) - V; 
bN = 0.125*exp(V/80);
end




function out = getpeaks(V)
[yc,xc] = findpeaks(V);
dim = size(yc);
  

coor = zeros(dim(1),2);
coor(:,2)=yc; 
coor(:,1)=xc;

out = [];

for i = 1:dim(1)
    if coor(i,2) >= -55
        out = [out;coor(i,:)];
    end
end


end


function sol = filt_peak(out,it)
sol = zeros(1,2);
sol(1,1) = it;
sol(1,2) = out(2,2);
disp(sol)
end


%%%%%% Brouillon formule capacité:

      %Cm(i+1) = Cm(i) * ((1 + 0.0048*(z(i+1) - z(i)))/(1 - 0.002*(z(i+1) - z(i)))); %Equation source Optocapacitance


      %Cm(i) = 0.824*1e-6 + (2.2*1e-6)/(31-z(i)); 
      %Cm(i) = (0.824*1e-6 + (2.2*1e-6)/(31-z(i)))*(d*l); 

      %Cm(i+1) = Cm(i)*(1 + 0.0048*(z(i+1)-z(i)))/(1 - 0.002*(z(i+1)-z(i))); 
      %Cm(i+1)= Cm(1)*(1 + 0.001*(z(i+1) - T_init)); 
