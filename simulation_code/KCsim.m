function [KC inh] = KCsim(KC,PN)

%compute excitatory input
KC.Vm       = KC.wPNKC' * kron(PN.rates,ones(PN.PN_per_glom,1));
KC.Vmspont  = KC.wPNKC' * kron(PN.spont,ones(PN.PN_per_glom,1));

%compute inhibitory input (assuming feedforward for now)
inh.odor    = KC.wInhKC * PN.wPNInh * PN.rates;
inh.spont   = KC.wInhKC * PN.wPNInh * PN.spont;

%divisive or subtractive inhibition?
if(KC.divisiveInh)
    KC.Vm        = KC.Vm./inh.odor;
    KC.Vmspont   = KC.Vmspont./inh.spont;
else
    KC.Vm        = KC.Vm - inh.odor;
    KC.Vmspont   = KC.Vmspont - inh.spont;
end


%compute the threshold
temp    = KC.Vm(:,1:110);
temp    = sort(temp(:),'descend');
KC.thr  = temp(length(temp) * KC.sparse);
KC.thr  = KC.thr * ones(KC.ncells,1);


%get active cells:
deltaI      = KC.Vm - KC.thr*ones(1,size(KC.Vm,2));
KC.rates    = KC.gain * (deltaI>0) .* tanh(deltaI/KC.gain);

%check for spontaneous activity:
deltaI      = KC.Vmspont - KC.thr;
KC.spont    = KC.gain*(deltaI>0).*tanh(deltaI/KC.gain);