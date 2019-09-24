realonly = 1;
ORN = {}; %olfactory receptor neurons
ORN         = buildORNs(ORN, realonly);
LN={};  %lateral neurons (inhbition between glomeruli)
LN          = buildLN(LN);
[PN, KC] = initialize_cells(ORN,realonly);

time = (PN.tmin+PN.dt:PN.dt:PN.tmax);
meanResp=[];clear PN_t;
for i=1:110
    disp(i);
    [~,PN_t,~,~] = getPNdynamics(ORN,PN,LN,i);
    meanResp(:,i) = mean(PN_t(:,time>PN.tOn & time<PN.tOff),2);
end
%%

figure(1);clf;
plot(bsxfun(@minus,meanResp,PN_t(:,1)),PN.rates_RW(:,1:110),'.')
hold on;
plot([0 200],[0 200],'r')
xlabel('dynamic model');
ylabel('RW model');
box off;
% an important thing to remember is that PN.rates is not the current model
% (I'm actually not too sure what it is... it's not a very good fit to the
% RW model.)