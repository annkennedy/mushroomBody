figure(1);clf

subplot(2,3,1);
image(corr(ORN.rates(ORN.HCList,1:110)')*32+32);axis equal;box off

subplot(2,3,2);
image(corr(PN.rates(ORN.HCList,1:110)')*32+32);axis equal;box off

subplot(2,3,3);
image(corr(PN.rates_RW(ORN.HCList,1:110)')*32+32);axis equal;box off

subplot(2,3,4)
image(corr(ORN.rates(ORN.HCList,1:110))*32+32);axis equal;box off

subplot(2,3,5)
image(corr(PN.rates(ORN.HCList,1:110))*32+32);axis equal;box off

subplot(2,3,6)
image(corr(PN.rates_RW(ORN.HCList,1:110))*32+32);axis equal;box off
%%
[C1,C2] = deal(zeros(110));
for rep = 1:10
    C = corr(experiment(1).results{rep}(:,1:110));
    C(isnan(C))=0;
    C = C-diag(diag(C))+eye(110);
    C1 = C1 + C/10;
    C = corr(experiment(2).results{rep}(:,1:110));
    C(isnan(C))=0;
    C = C-diag(diag(C))+eye(110);
    C2 = C2 + C/10;
end

figure(1);clf;
subplot(1,2,1);
image(C1*32+32);axis equal;axis tight;axis off
hold on;

for i = unique(ORN.odorclass(1:110))
    plot(get(gca,'xlim'),.5+[1 1]*sum(ORN.odorclass<i),'k');
    plot(.5+[1 1]*sum(ORN.odorclass<i),get(gca,'ylim'),'k');
end

subplot(1,2,2);
image(C2*32+32);axis equal;axis tight;axis off
hold on;

for i = unique(ORN.odorclass(1:110))
    plot(get(gca,'xlim'),.5+[1 1]*sum(ORN.odorclass<i),'k');
    plot(.5+[1 1]*sum(ORN.odorclass<i),get(gca,'ylim'),'k');
end

%%
[rates1,rates2] = deal([]);
for rep = 1:10
    rates1 = [rates1; sum(experiment(1).results{rep}(:,1:110)~=0)];
    rates2 = [rates2; sum(experiment(2).results{rep}(:,1:110)~=0)];
end

figure(1);clf;
subplot(1,2,1);hold on;
bar(mean(rates1));
errorbar(mean(rates1),std(rates1),'k.');
ylim([0 1100]);

subplot(1,2,2);hold on;
bar(mean(rates2));
errorbar(mean(rates2),std(rates2),'k.');
ylim([0 1100]);
%%

bins = linspace(0,1100,20);
n1=hist(rates1(:),bins);
n2=hist(rates2(:),bins);

figure(2);clf;hold on;
area(bins,n2,'facealpha',.5)
area(bins,n1,'facealpha',.5)