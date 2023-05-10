function out = test_model_T(Temp_max, temps_monte, t_exp, start_time, Temp_init)
    
    temps_monte = temps_monte * 100;
    start_time = start_time * 100;
    pente = Temp_max/temps_monte;
    tau = 0.1* 1e2;

    y_monte = @(t) pente*t + Temp_init - pente*start_time;
    y_descente = @(t) Temp_max*(exp(-(t-temps_monte-start_time)/tau)) + Temp_init;

    out = [];
    out = [out, 0];

    for i=2:length(t_exp)
        if i <= start_time
            out = [out, Temp_init];
        end

        if  i > start_time && i <= (temps_monte+start_time) 
            value = y_monte(i); 
            out = [out, value];
        end

        if i > (temps_monte+start_time) 
            value = y_descente(i); 
            out = [out, value];
        end


        
    end

    plot(out)
end

