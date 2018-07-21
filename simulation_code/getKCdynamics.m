function [KC_s KC_t inhstore] = getKCdynamics(PN,nreps,KC)

dt      = KC.dt;
time    = KC.tmin+dt:dt:KC.tmax;
tau_inh = 0.05; %assume APL inhibition is slow (following the lead of LN inhibition)
tau_Is  = 0.01;% synaptic time constant, KC->inh

thr = KC.thr + KC.spont;

w = KC.wPNKC;

for rep = 1:nreps
    Vm = zeros(KC.ncells,length(time));
    sp = zeros(KC.ncells,length(time));
    PN_s = smoothts([zeros(size(PN{1},1),1) PN{rep}],'e',KC.tau_s/KC.dt); %convolve with synaptic filter
    inh = zeros(1,length(time));
    Is  = zeros(1,length(time));
    for t = 2:length(time)
        dIsdt   = -Is(t-1) + KC.wKCInh*sp(:,t-1)*1e4;
        dinhdt  = -inh(t-1) + Is(t-1);
        dKCdt   = -Vm(:,t-1) + w'*PN_s(:,t) - (KC.wInhKC*inh(t-1));
        
        Vm(:,t) = Vm(:,t-1) + dKCdt*dt/KC.taum;
        inh(t)  = inh(t-1) + dinhdt*dt/tau_inh;
        Is(t)   = Is(t-1) + dIsdt*dt/tau_Is;

        sp(Vm(:,t)>thr,t) = 1;
        Vm(Vm(:,t)>thr,t) = 0;

    end
    KC_s{rep} = sparse(sp);
    KC_t{rep} = Vm;
    inhstore(rep,:) = inh;
end