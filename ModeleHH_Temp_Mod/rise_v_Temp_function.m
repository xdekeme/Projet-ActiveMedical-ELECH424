
function v_list = rise_v_Temp_function(Tm_init)
    %Define global constants
    v_l = 10.589; %in mV
    V_rest=-70; %in mv

    v_list = [V_rest];

    gl = 0.3; %in /kohm*cm^2 **also gmax_l
    %c0 = 0.824; %in microfarad/cm^2
    dt = 0.001; %in ms
    gna=120.0;%in /kohm*cm^2 **also gmax_l
    gk=36.0;%in /kohm*cm^2 **also gmax_l
    v_na=115.0 ;%in mV
    v_k=-12.0;%in mV
    d=2;
    d=d/10000; %in cm **axon_diameter
    dx=0.025; %in cm              on a 36 points répartis sur 0.9 cm et chacun séparé de 0.025 cm
    ro_i = .0345; %in KOhm*cm
    %ro_i = .0145; %in KOhm*cm
    ro_e = .3; %in KOhm*cm
    x0 = 0.1; %in cm
    z0 = 0.1; %in cm
    x1 = 0.6; %in cm
    z1 = 0.1; %in cm
    global R;
    R=(4*ro_i*dx*dx)/d;
    global T;
    T=25;

    T_list = [T];

    %k1=3.0^((T-6.3)/10);

    %create matricies bounded by time and space inputs
    %npoints = input('length of the axon (cm): ');
    npoints=0.9;
    npoints=npoints/dx;
    % stop= input('period of time to run test(ms): '); %%%%%%
    stop = 16;
    %stop=0.125;
    stop=stop/dt;
    M=zeros(npoints+1,stop);
    H=zeros(npoints+1,stop);
    N=zeros(npoints+1,stop);
    kv1=zeros(npoints+2,stop);
    kv2=zeros(npoints+2,stop);
    kv3=zeros(npoints+2,stop);
    kv4=zeros(npoints+2,stop);
    km1=zeros(npoints+1,stop);
    km2=zeros(npoints+1,stop);
    km3=zeros(npoints+1,stop);
    km4=zeros(npoints+1,stop);
    kn1=zeros(npoints+1,stop);
    kn2=zeros(npoints+1,stop);
    kn3=zeros(npoints+1,stop);
    kn4=zeros(npoints+1,stop);
    kh1=zeros(npoints+1,stop);
    kh2=zeros(npoints+1,stop);
    kh3=zeros(npoints+1,stop);
    kh4=zeros(npoints+1,stop);
    V=zeros(npoints+2,stop);
    Ve=zeros(npoints+2,stop);
    TT=zeros(npoints+2,stop+1);
    TT1=zeros(npoints+2,stop+1);
    DTT=zeros(npoints+2,stop);
    K1=zeros(npoints+2,stop);
    CC=zeros(npoints+2,stop);
    dCC=zeros(npoints+2,stop);
    Ic=zeros(npoints+2,stop);

    %test stimulation current
    I= zeros(1,stop);
    %Amp= input('input test stimulation intensity(uA):');
    Amp=0;
    Pw=0.1; %pulsewidth in ms
    Pw=Pw/dt;
    St=0.1/dt;
    for i=St:St+Pw
        I(i)=-1*Amp;
    end

    %block stimulation current
    Ib= zeros(1,stop);
    %Ampb= input('input block stimulation intensity(uA):');

    %Frqb= input('input block stimulation frequency(kHz):');

    for i=1:stop
        Ib(i)=0;
    end
    %Prdb=1/(Frqb*dt);
    %Prdb=round(Prdb/2);
    %for i=(Prdb+1):Prdb*2:stop
    %   for j=i:(i+Prdb-1)
    %    Ib(j)=Ampb;
    %    end
    %end

    %calculate external voltages
    ts=1;
    while ts<=stop,
        for step=1:npoints+2,
            Ve(step,ts)=(ro_e*I(ts))/(4*pi*((((step-2)*dx-x0)^2)+z0^2)^0.5);
            Ve(step,ts)= Ve(step,ts)+(ro_e*Ib(ts))/(4*pi*((((step-2)*dx-x1)^2)+z1^2)^0.5);
            
        end
        ts=ts+1;
    end

    %define initial conditions at time 1 for v and gate variables
    V(:,1)=0;
    M(:,1)=0.053;
    H(:,1)=0.596;
    N(:,1)=0.318;
    TT(:,:)= T;
    %Tr = input('rising time for temperature stim(ms): '); %%%%%%
    Tr = 1;
    Tr = Tr/dt;
    %Tw = input('width for temperature stim(cm): '); %%%%%%%
    Tw = 0.05;
    Tw = Tw/dx;
    %Tm = input('maximum temperature change(C): '); %%%%%
    Tm = Tm_init;
    disp(Tm)
    %k = input('temperature coefficiency(KuF/cm^2): '); %%%%%
    k = 2.2;
    Tc=31;
    c0=1-k/(Tc-18.5);

    ts=20;

    disp("1")
    % Code to create the temperature vector
    while ts<=stop,
        for step=1:npoints+1,
            TT1(step,ts)=Tm_init*exp(-(step-0.45/dx)^2/(2*Tw^2));   % fonction qui code pour la diffusion de la température dans le matériau
        end


        ts=ts+1;
    end


    disp("2")
    ts=20;
    while ts<=(Tr+20),
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*(ts-20)/(Tr);  % en un point la variation temporelle de T est linéaire
        end

        T_list = [T_list, TT(round(npoints/2),ts)];

        ts=ts+1;
    end
    while ts<=stop,
        for step=1:npoints+1,
            TT(step,ts)=T+TT1(step,ts)*exp(-(ts-(Tr+20))/(100/dt));
        end

        T_list = [T_list, TT(round(npoints/2),ts)];

        ts=ts+1;
    end

    for step=1:npoints+1,
        TT(step,stop+1)= TT(step,stop);
    end

    disp("3")
    %create the sensitivity to temperature
    for i=2:npoints+1
        for j=2:stop+1
            K1(i-1,j-1)=3^((TT(i-1,j-1)-6.3)/10);
            DTT(i-1,j-1)=(TT(i-1,j)-TT(i-1,j-1))/dt;  % dérivée temporelle de la température
            CC(i-1,j-1)=c0+(k/(Tc-TT(i-1,j-1)));
            dCC(i-1,j-1)=k*DTT(i-1,j-1)/((Tc-TT(i-1,j-1))^2); % dérive temporelle de la capa
        end
    end

    t=2;
    disp("4")
    while t<=stop,
    	%start caculating k1
    	for step=2:npoints+1;
            v=V(step,t-1);
            v1=V(step-1,t-1);
            v3=V(step+1,t-1);
            ve=Ve(step,t-1);
            ve1=Ve(step-1,t-1);
            ve3=Ve(step+1,t-1);
            m=M(step-1,t-1);
            n=N(step-1,t-1);
            h=H(step-1,t-1);
            k1=K1(step-1,t-1);
            C=CC(step-1,t-1);
            dC=dCC(step-1,t-1);
            [kv,I_na,I_k,I_l,Ii]=fv1(v1,v,v3,ve1,ve,ve3,m,n,h,R,C,dC);

            kv1(step-1,t-1)=kv;
            Ina(step-1,t-1)=I_na;
            Ik(step-1,t-1)=I_k;
            Il(step-1,t-1)=I_l;
            I_in(step-1,t-1)=Ii;
            Ic(step-1,t-1)=CC(step-1,t-1)*((V(step-1,t)-V(step-1,t-1))/dt)+...
                (V(step-1,t-1)-70)*dCC(step-1,t-1);
            km1(step-1,t-1)=fm(v,m,k1);
            kn1(step-1,t-1)=fn(v,n,k1);
            kh1(step-1,t-1)=fh(v,h,k1);
        end

        %start caculate k2
        for step=2:npoints+1;
            k1=K1(step-1,t-1);
            C=CC(step-1,t-1);
            dC=dCC(step-1,t-1);
            v=V(step,t-1)+dt*kv1(step,t-1)/2;
            v1=V(step-1,t-1)+dt*kv1(step-1,t-1)/2;
            v3=V(step+1,t-1)+dt*kv1(step+1,t-1)/2;
            ve=Ve(step,t-1);
            ve1=Ve(step-1,t-1);
            ve3=Ve(step+1,t-1);
            m=M(step-1,t-1)+dt*km1(step-1,t-1)/2;
            n=N(step-1,t-1)+dt*kn1(step-1,t-1)/2;
            h=H(step-1,t-1)+dt*kh1(step-1,t-1)/2;
            kv2(step-1,t-1)=fv(v1,v,v3,ve1,ve,ve3,m,n,h,R,C,dC);
            km2(step-1,t-1)=fm(v,m,k1);
            kn2(step-1,t-1)=fn(v,n,k1);
            kh2(step-1,t-1)=fh(v,h,k1);
        end
        %start caculate k3
        for step=2:npoints+1;
            k1=K1(step-1,t-1);
            C=CC(step-1,t-1);
            dC=dCC(step-1,t-1);
            v=V(step,t-1)+dt*kv2(step,t-1)/2;
            v1=V(step-1,t-1)+dt*kv2(step-1,t-1)/2;
            v3=V(step+1,t-1)+dt*kv2(step+1,t-1)/2;
            ve=Ve(step,t-1);
            ve1=Ve(step-1,t-1);
            ve3=Ve(step+1,t-1);
            m=M(step-1,t-1)+dt*km2(step-1,t-1)/2;
            n=N(step-1,t-1)+dt*kn2(step-1,t-1)/2;
            h=H(step-1,t-1)+dt*kh2(step-1,t-1)/2;
            kv3(step-1,t-1)=fv(v1,v,v3,ve1,ve,ve3,m,n,h,R,C,dC);
            km3(step-1,t-1)=fm(v,m,k1);
            kn3(step-1,t-1)=fn(v,n,k1);
            kh3(step-1,t-1)=fh(v,h,k1);
        end
        %start caculate k4
        for step=2:npoints+1;
            k1=K1(step-1,t-1);
            C=CC(step-1,t-1);
            dC=dCC(step-1,t-1);
            v=V(step,t-1)+dt*kv3(step,t-1);
            v1=V(step-1,t-1)+dt*kv3(step-1,t-1);
            v3=V(step+1,t-1)+dt*kv3(step+1,t-1);
            ve=Ve(step,t-1);
            ve1=Ve(step-1,t-1);
            ve3=Ve(step+1,t-1);
            m=M(step-1,t-1)+dt*km3(step-1,t-1);
            n=N(step-1,t-1)+dt*kn3(step-1,t-1);
            h=H(step-1,t-1)+dt*kh3(step-1,t-1);
            kv4(step-1,t-1)=fv(v1,v,v3,ve1,ve,ve3,m,n,h,R,C,dC);
            km4(step-1,t-1)=fm(v,m,k1);
            kn4(step-1,t-1)=fn(v,n,k1);
            kh4(step-1,t-1)=fh(v,h,k1);
        end
        %caculate v,m,n,h,p
        for step=1:npoints;
            V(step+1,t)=V(step+1,t-1)+dt*(kv1(step,t-1)+2*kv2(step,t-1)+...
                2*kv3(step,t-1)+kv4(step,t-1))/6;
            M(step,t)=M(step,t-1)+dt*(km1(step,t-1)+2*km2(step,t-1)+...
                2*km3(step,t-1)+km4(step,t-1))/6;
            N(step,t)=N(step,t-1)+dt*(kn1(step,t-1)+2*kn2(step,t-1)+...
                2*kn3(step,t-1)+kn4(step,t-1))/6;
            H(step,t)=H(step,t-1)+dt*(kh1(step,t-1)+2*kh2(step,t-1)+...
                2*kh3(step,t-1)+kh4(step,t-1))/6;


            if step == round(npoints/2);
                v_list = [v_list, V(2,t)+V_rest];
                %disp(v_list)
            end

        end
        V(1,t)=V(2,t);V(npoints+2,t)=V(npoints+1,t);
        %   plot(V(:,t))
        %     axis([1,npoints+2,-100,120]);
        %     drawnow;
        t=t+1;
    end
    %end of updating v, m, and h
    %plot(v(:,:))

    section_long=ones(1,stop);




   % plot of temp at a certain time
    % figure;
    % plot(TT(:,3500))
  % 
  % %   plot of temp at a certain place
  %   figure;
  %   plot(T_list);
    % figure
    % plot(CC(npoints/2,:))
 
end

