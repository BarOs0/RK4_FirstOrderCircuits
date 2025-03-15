function y = RK4(t, y0, R, C, L, h, A, F, option)

    N = length(t);
    y = zeros(1,N);
    y(1) = y0;

    for n = 1:N-1
        k1 = h.*Fout(t(n), y(n), R, C, L, A, F, option);
        k2 = h.*Fout(t(n) + h/2, y(n) + k1/2, R, C, L, A, F, option);
        k3 = h.*Fout(t(n) + h/2, y(n) + k2/2, R, C, L, A, F, option);
        k4 = h.*Fout(t(n) + h, y(n) + k3, R, C, L, A, F, option);
        dyn = (1/6).*(k1 + 2.*(k2 + k3) + k4);
        y(n+1) =  y(n) + dyn;
    end
end

