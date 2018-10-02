
[~,i] = sort(ORN.spont(ORN.HCList),'ascend');

figure(1);clf;
subplot(3,1,1);
nCxns = hist(gList,1:55);
bar(nCxns(ORN.HCList(i))/sum(nCxns(ORN.HCList)))
title('proportion of PN-KC cxns');
box off;

subplot(3,1,2);
bar(ORN.spont(ORN.HCList(i)));
title('spont ORN firing rate (Hz)')
box off;

subplot(3,1,3);
bar(mean(ORN.rates(ORN.HCList(i),1:110),2))
hold on;
errorbar(1:23, mean(ORN.rates(ORN.HCList(i),1:110),2),...
               std(ORN.rates(ORN.HCList(i),1:110),[],2)/sqrt(110),'k.')
title('evoked ORN firing rate (Hz)');
box off;


%%
spont   = ORN.spont(ORN.HCList(i));
cxns    = nCxns(ORN.HCList(i))'/sum(nCxns(ORN.HCList));
evoked  = mean(ORN.rates(ORN.HCList(i),1:110),2);

%%
figure(2);clf
subplot(1,2,1);

mdl = fitlm(table(spont,cxns));
h = plot(mdl);
set(h(1),'marker','.','markersize',15,'color',[0 0 0]);
box off;


subplot(1,2,2)
mdl = fitlm(table(evoked,cxns));
h = plot(mdl);
set(h(1),'marker','.','markersize',15,'color',[0 0 0]);
box off



