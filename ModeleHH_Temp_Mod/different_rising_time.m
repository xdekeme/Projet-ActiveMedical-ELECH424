
T=18.5;

dx=0.025;
npoints=0.9;
npoints=npoints/dx;

stop =7;
dt = 0.001;
stop=stop/dt;

TT=zeros(npoints+2,stop+1);
TT1=zeros(npoints+2,stop+1);

Tm = 8;

Tr_list = [0.01 0.05 0.1 0.5 1.0 2.0];

Tw = 0.05;
Tw = Tw/dx;

for i = 1:length(Tr_list)

    T_list = [T];

    ts=20;

    while ts<=stop,
      for step=1:npoints+1,
          TT1(step,ts)=Tm*exp(-(step-0.45/dx)^2/(2*Tw^2));
      end
      ts=ts+1;
    end

    ts=20;
    while ts<=((Tr_list(i)/dt)+20),
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*(ts-20)/((Tr_list(i)/dt));
        end

    

    T_list = [T_list, TT(round(npoints/2),ts)];

    ts=ts+1;
    end

    while ts<=stop,
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*exp(-(ts-((Tr_list(i)/dt)+20))/(100/dt));
        end

    

    T_list = [T_list, TT(round(npoints/2),ts)];

    ts=ts+1;
    end

    for step=1:npoints+1,
        TT(step,stop+1)= TT(step,stop);
    end

    if i == 1
       figure;
       plot(T_list); 
       hold on 
    end

    plot(T_list);
    hold on

end








