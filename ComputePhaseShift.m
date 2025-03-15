function ps = ComputePhaseShift(signal1, signal2, F, t)

    sgn1 = signal1(floor(end/2):end);
    sgn2 = signal2(floor(end/2):end);

    zero1 = find(sgn1(1:end-1) < 0 & sgn1(2:end) >= 0, 1);
    zero2_candidates = find(sgn2(1:end-1) < 0 & sgn2(2:end) >= 0);
    zero2 = zero2_candidates(zero2_candidates > zero1);

    if ~(isempty(zero2))
        zero2 = zero2(1);
    else
        ps = 90; return;
    end

    dt = abs(t(zero2) - t(zero1));
    ps = 360 - (360 * F * dt);

    if ps >= 180
        ps = ps - 360;
    end

end