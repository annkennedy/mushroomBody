function PN = PNsim(PN,ORN,realonly)

nodors = size(ORN.rates,2);

%we don't know what spontaneous PN activity looks like, so I'll just assume
%it takes the form of the ORN spont rates scaled by the lateral inhibition.
PN.spont    = ORN.spont * PN.inhsc/(sum(ORN.spont) + PN.inhadd);

for i = 1:nodors

    if(realonly)
        b = sum(ORN.rates(:,i) + ORN.spont)*PN.ncells/length(PN.HCList);
    else
        b = sum(ORN.rates(:,i) + ORN.spont);
    end
    inh = PN.inhsc/(b + PN.inhadd);

    PN.rates(:,i)   = PN.Rmax_dyn * tanh((ORN.rates(:,i)+PN.offset)*PN.tanhsc*inh/PN.Rmax_dyn) + PN.spont; %the new PN model
    PN.rates(:,i)   = max(PN.rates(:,i),0);

%adding in the old model to test stuff
    b_old = sum(ORN.rates(:,i))*PN.ncells/length(PN.HCList);
    PN.rates_RW(:,i)   = (165*ORN.rates(:,i).^1.5)./(12.^1.5 + ORN.rates(:,i).^1.5 + (.05*b_old).^1.5);
    PN.rates_RW(:,i)   = real(PN.rates_RW(:,i).*(ORN.rates(:,i)>=0)); %zeros PNs with negative input



end

if(realonly)
    PN.rates(setdiff(1:PN.ncells,ORN.HCList),:)   = 0;
    PN.spont(setdiff(1:PN.ncells,ORN.HCList))     = 0;
end 