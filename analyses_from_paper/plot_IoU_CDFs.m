

h=[];legStrs={};
figure(1);clf;
% figure(2);clf;
% figure(3);clf;
for i=3%:10%:length(experiment)
    opts = experiment(i);
    OLstack = [];
    silentBars = [];
silentCellBars = [];
    for rep = 1:opts.nreps
        M = opts.results{rep}(:,1:110);
        Mand    = double(M~=0)'*double(M~=0);
%         Mor     = sum(M~=0)'+sum(M~=0) - Mand;
        Mor = ones(110,1)*sum(M~=0);
        overlap = Mand./Mor;% + tril(nan(size(Mand)));
        overlap = overlap(~isnan(overlap))';
        OLstack = [OLstack overlap];
        silentBars(end+1) = sum(sum(M~=0)==0);
        silentCellBars(end+1) = sum(sum(M'~=0)==0);
    end

    figure(1);hold on;
    h(end+1) = plotCDFFit(overlap);
    figure(2);subplot(2,1,1);hold on;
    bar(i,mean(silentBars))
    errorbar(i,mean(silentBars),std(silentBars),'k');
    subplot(2,1,2);hold on;
    bar(i,mean(silentCellBars))
    errorbar(i,mean(silentCellBars),std(silentCellBars),'k');
    
    
    str = ['sparsity=' num2str(opts.sparsity) ' | '];
    if(opts.realonly),    str = [str 'real-only | ']; end
    if(opts.useAPL),      str = [str '+APL | ']; else, str = [str 'no APL | ']; end
    if(opts.homeostatic), str = [str 'homeo on | ']; end
    if(opts.shuffleORNs), str = [str 'shuffle ORNs | ']; end
    if(opts.uniformPNKC), str = [str 'uniform wPNKC | ']; end
    if(opts.useStaticPN), str = [str 'RW PN model']; end
    legStrs{end+1} = str;
                                        
end

figure(1);
%%
% 10% of cells for each odor
for i=1:10000
    vals = rand(1,2000);
    s = sort(vals);
    M(:,i) = vals>s(round(2000*.9));
end
Mand    = double(M~=0)'*double(M~=0);
% Mor     = sum(M~=0)'+sum(M~=0) - Mand;
Mor = ones(10000,1)*sum(M~=0);
overlap = Mand./Mor;% + tril(nan(size(Mand)));
overlap = overlap(~isnan(overlap))';
%%
h(end+1)=plotCDFFit(overlap);
legStrs{end+1} = 'rand vectors';
set(h(end),'linewidth',2,'color','k');


xlabel('fraction of coactive neurons');
ylabel('fraction of odor pairs');
% set(gca,'xscale','log')
legend(h,legStrs);
xlim([1e-3 1])

%%
figure(2);
subplot(2,1,1);
title('silent odors');
subplot(2,1,2);
title('silent cells');
