function H = Transmitance(R,L,C,freq,option)

   w = 2*pi*freq;
   s = 1j*w;

    switch (option)

        case('RC')
            H = 1./(R*C*s + 1);
        case('CR')
            H = (R*C*s)./(R*C*s + 1);
        case('RL')
            H = (L*s)./(L*s + R);
        case('LR')
            H = R./(L*s + R);
    end
end

