function vpp = ComputeVpp(signal)
    
    signal = signal(floor(end/2):end); % established signal values
    dc_removed = signal - mean(signal); % dc offset removed
    vpp = abs(max(dc_removed)-min(dc_removed)); % calculating vpp
end

