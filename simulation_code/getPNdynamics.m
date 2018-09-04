function [ORN_t,PN_t,LNR,PN_RW] = getPNdynamics(ORN,PN,LN,stimodor)

dt      = PN.dt;
time    = PN.tmin-1+dt:dt:PN.tmax; % we simulate for an extra bit at the beginning to get to steady state
stim    = (time>PN.tOn)&(time<PN.tOff);
odor    = smoothts(ORN.spont*ones(size(stim)) + ORN.rates(:,stimodor)*stim,'e',.020/dt);   % odor input to ORNs


ORN_t   = ORN.spont * ones(size(stim));                                                    % ORN timecourse
PN_t    = PN.spont * ones(size(time));
LN_t    = 300*ones(size(stim));
LNR     = ones(size(stim));
inhA    = 50*ones(size(stim));
inhB    = 50*ones(size(stim));
inh_PN = 0;
inh_LN = 0;

PN_RW    = PN.spont * ones(size(time));
% PN_inhKO = PN.spont*ones(size(time));

for t = 2:length(time)
    ORN_delta = ORN_t(:,t-1) - ORN.spont; %ORNs relative to spontaneous (use for PN tanh nonlinearity) (drop this?)

    dinhAdt = -inhA(t-1) + LNR(t-1);
    dinhBdt = -inhB(t-1) + LNR(t-1);

    dORNdt = -ORN_t(:,t-1) + odor(:,t);
    dPNdt  = -PN_t(:,t-1)  + PN.spont + 200*tanh((ORN_delta + PN.offset)*PN.tanhsc/200*inh_PN);
    dLNdt  = -LN_t(t-1)    + mean(ORN_t(:,t-1))^(3)*51/23/2*inh_LN;



    inhA(t)= inhA(t-1) + dinhAdt*dt/LN.tauGA;
    inhB(t)= inhB(t-1) + dinhBdt*dt/LN.tauGB;

    inh_PN = PN.inhsc/(PN.inhadd + .25*inhA(t) + .75*inhB(t));
    inh_LN = LN.inhsc/(LN.inhadd + inhA(t));

    LN_t(t)     = LN_t(t-1) + dLNdt*dt/LN.taum;
    LNR(t)      = (LN_t(t) - LN.thr)*(LN_t(t)>LN.thr);
    ORN_t(:,t)  = ORN_t(:,t-1) + dORNdt*dt/ORN.taum;
    PN_t(:,t)   = PN_t(:,t-1) + dPNdt*dt/PN.taum;
    PN_t(:,t)   = max(PN_t(:,t),0);
    
    inh_static = PN.inhsc/(PN.inhadd + sum(ORN_t(:,t-1)));
    dPNRWdt    = -PN_RW(:,t-1)    + PN.spont + max(200*tanh((ORN_delta+PN.offset)*PN.tanhsc/200*inh_static),0);
    PN_RW(:,t) = PN_RW(:,t-1) + dPNRWdt*dt/PN.taum;
    
%     inhKO      = PN.inhsc/(PN.inhadd+sum(ORN.spont));
%     dPNinhKOdt = -PN_inhKO(:,t-1) + PN.spont + max(200*tanh((ORN_delta+PN.offset)*PN.tanhsc/200*inhKO),0);
%     PN_inhKO(:,t) = PN_inhKO(:,t-1) + dPNinhKOdt*dt/PN.taum;

end

nonHC = setdiff(1:51,ORN.HCList);
if(sum(ORN.spont(nonHC))==0) %a cheap hack to indirectly test the value of realonly
    for i=1:length(nonHC)
        PN_t(nonHC(i),:)=0; %get rid of cells we aren't simulating- background effects of inh could muck up future analysis
        PN_RW(nonHC(i),:)=0;
%         PN_inhKO(nonHC(i),:)=0;
    end
end

% cut out the first second of simtime (used to get the model to steady-state spontaneous conditions)
ORN_t = ORN_t(:,1/dt+1:end);
PN_t = PN_t(:,1/dt+1:end);
PN_RW = PN_RW(:,1/dt+1:end);
% PN_inhKO = PN_inhKO(:,1/dt+1:end);