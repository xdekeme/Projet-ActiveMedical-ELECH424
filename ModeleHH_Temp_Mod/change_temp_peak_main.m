T=25;



dx=0.025;
npoints=0.9;
npoints=npoints/dx;
v_list_tot = [];

stop =7;
dt = 0.001;
stop=stop/dt;

TT=zeros(npoints+2,stop+1);
TT1=zeros(npoints+2,stop+1);

%Tm_list = [7 7.25 7.5 7.75 8];
%Tm_list = [7.95 7.96 7.97 7.98 7.99]; %The model starts acting at 8°C jump
%for a rising time of 1 and a starting temp of 18.5°C
Tm_list = [4 4.5 5 5.5]; 

Tr = 1;
%Tr = Tr/dt;

Tw = 0.05;
Tw = Tw/dx;


for i = 1:length(Tm_list)

    T_list = [T];
    ts=20;

    while ts<=stop,
      for step=1:npoints+1,
          TT1(step,ts)=Tm_list(i)*exp(-(step-0.45/dx)^2/(2*Tw^2));
      end
      ts=ts+1;
    end

    ts=20;
    while ts<=((Tr/dt)+20),
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*(ts-20)/((Tr/dt));
        end

    T_list = [T_list, TT(round(npoints/2),ts)];

    ts=ts+1;
    end

    while ts<=stop,
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*exp(-(ts-((Tr/dt)+20))/(100/dt));
        end

    T_list = [T_list, TT(round(npoints/2),ts)];

    ts=ts+1;
    end

    for step=1:npoints+1,
        TT(step,stop+1)= TT(step,stop);
    end
    
    % figure;
    % plot(T_list);
    % title("Evolution of temperature at different temperature jump (at starting temperature: " + T + "°C)")
    % ylabel("Temperature (°C)")
    % xlabel("Time (us)")
    % legend("1°C", "3°C", "5°C", "7°C", "9°C", "11°C")
    % hold on


end

figure;
hold on;

for i = 1:length(Tm_list)
    v_list = rise_v_Temp_function(Tm_list(i));
    plot(v_list);
end

title("AP at different temperature jump (at starting temperature: " + T + "°C)")
ylabel("Voltage (mV)")
xlabel("Time (us)")
legend(Tm_list(1) + "°C", Tm_list(2) +"°C", Tm_list(3) +"°C", Tm_list(4) + "°C");

figure;
plot(capa_Temp_function(Tm_list(1)));
xlabel("Time (us)");
ylabel("Capacitance membrane (uF/cm^{2})");
title("Value of the capacitance for an initial temperature of the axon 25°C and a temperature peak of 10°C")

% figure;
% plot(capa_Temp(Tm_list(1)));
% xlabel("Time (us)");
% ylabel("Temperature (°C)");
% title("Different temperature peaks in a row")
    