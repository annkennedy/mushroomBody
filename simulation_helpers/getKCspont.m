function spont = getKCspont(PN, nreps, KC)

dt      = KC.dt;
time    = KC.tmin+dt:dt:KC.tOn; %hardcodes in stim onset at 0, not the best

w = KC.wPNKC;

Vm_mean = zeros(KC.ncells,length(time));
for rep = 1:nreps
    Vm = zeros(KC.ncells,length(time));
    PN_s = smoothts([zeros(size(PN{1},1),1) PN{rep}],'e',KC.tau_s/dt);
    for t = 2:length(time)
        dKCdt       = -Vm(:,t-1) + w'*PN_s(:,t);
        Vm(:,t)   = Vm(:,t-1) + dKCdt*dt/KC.tau_m;
    end
    Vm_mean = Vm_mean + Vm/nreps;
end

spont = mean(Vm_mean(:,250:end),2);