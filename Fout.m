function ddt = Fout(t, y, R, C, L, A, F, option)

% y - Voltage or current, depends of case

Vin = Fnk(A,F,t);

    switch(option)
        case('c')
            ddt = (Vin - y)/(R*C);
        case('l')
            ddt = (Vin - R*y)/L;
    end
end

