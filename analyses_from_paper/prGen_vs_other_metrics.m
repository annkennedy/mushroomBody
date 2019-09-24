

% generalizationStats(generalizationStats==0)=1;
% prGen = squeeze(mean(generalizationStats,1));


opts = experiment(1);
prGen = prGen_full;

overlap = zeros(110);
for rep = 1:opts.nreps
    M = opts.results{rep}(:,1:110);
    Mand    = double(M~=0)'*double(M~=0);
    Mtest = ones(110,1)*sum(M~=0);
    overlap = overlap + (Mand./Mtest)/opts.nreps;
%     Mor     = sum(M~=0)'+sum(M~=0) - Mand;
%     overlap = overlap + (Mand./Mor)/opts.nreps;
end
overlap = Cpn;

clear meanOL;
OLrange = -0.5:.025:1;
for i=1:length(OLrange)-1
    inds = (overlap>=OLrange(i)) & (overlap<OLrange(i+1));
    meanOL(i) = 1-(nanmean(prGen(inds))+1)/2;
    stdOL(i)  = nanstd(prGen(inds)/2);
end

figure(1);clf;
plot(overlap,1-(prGen+1)/2,'b.');
hold on;
drawvar(OLrange(1:end-1)+.025/2,meanOL,'k',1,stdOL)
ylabel('probability of overgeneralization');
xlabel('overlap with target odor');
xlim([-0.5 1]);box off;
ylim
set(gca,'xtick',-0.5:.25:1);
%%

nSamp = 10;
[broad, narrow] = find_public_and_private_odors(ORN,nSamp);

Cpn         = corr(PN.rates(ORN.HCList,1:110)) - eye(110);
useX = Cpn;

figure(2);clf;
subplot(1,2,1);hold on;
prGen = prGen_full;

plot(useX,1-(prGen+1)/2,'.','color',[.75 .75 .75]);
plot(useX(broad,broad),1-(prGen(broad,broad)+1)/2,'o','markersize',5,'markeredgecolor','none','markerfacecolor','r')
plot(useX(narrow,narrow),1-(prGen(narrow,narrow)+1)/2,'o','markersize',5,'markeredgecolor','none','markerfacecolor','b')

subplot(1,2,2);hold on;
prGen = prGen_homeo;

plot(useX,1-(prGen+1)/2,'.','color',[.75 .75 .75]);
plot(useX(broad,broad),1-(prGen(broad,broad)+1)/2,'o','markersize',5,'markeredgecolor','none','markerfacecolor','r')
plot(useX(narrow,narrow),1-(prGen(narrow,narrow)+1)/2,'o','markersize',5,'markeredgecolor','none','markerfacecolor','b')
%%
class_Gen = zeros(10);
for i=1:10
    inds = ORN.odorclass==i;
    for j=1:10
        indsj = ORN.odorclass==j;
        class_Gen(i,j) = mean(mean(prGen(inds,indsj)));
    end
end
figure(1);clf;
image((1-(class_Gen+1)/2)*64*5);
colorbar
box off;
xlabel('test class');
ylabel('train class');








