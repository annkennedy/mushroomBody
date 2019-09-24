
nSamp = 10;
[broad, narrow] = find_public_and_private_odors(ORN,nSamp);

figure(1);clf;
for useExpt = 1:2
    Cpn         = corr(PN.rates(ORN.HCList,1:110)) - eye(110);
    Cpn_broad   = corr(PN.rates(ORN.HCList,broad)) - eye(nSamp);
    Cpn_narrow  = corr(PN.rates(ORN.HCList,narrow)) - eye(nSamp);

    [Ckc,overlap] = deal(zeros(size(Cpn)));
    [Ckc_broad,Ckc_narrow,overlap_broad,overlap_narrow] = deal(zeros(nSamp));
    for rep = 1:experiment(useExpt).nreps
        Ckc         = Ckc        + corr(experiment(useExpt).results{rep}(:,1:110))/experiment(useExpt).nreps;
        Ckc_broad   = Ckc_broad  + corr(experiment(useExpt).results{rep}(:,broad))/experiment(useExpt).nreps;
        Ckc_narrow  = Ckc_narrow + corr(experiment(useExpt).results{rep}(:,narrow))/experiment(useExpt).nreps;

        M       = experiment(useExpt).results{rep}(:,1:110);
        Mand    = double(M~=0)'*double(M~=0);
        Mor = ones(110,1)*sum(M~=0);
%         Mor     = sum(M~=0)'+sum(M~=0) - Mand;
        overlap = overlap + Mand./Mor/experiment(useExpt).nreps;

        M       = experiment(useExpt).results{rep}(:,broad);
        Mand    = double(M~=0)'*double(M~=0);
        Mor = ones(nSamp,1)*sum(M~=0);
%         Mor     = sum(M~=0)'+sum(M~=0) - Mand;
        overlap_broad = overlap_broad + Mand./Mor/experiment(useExpt).nreps;

        M       = experiment(useExpt).results{rep}(:,narrow);
        Mand    = double(M~=0)'*double(M~=0);
        Mor = ones(nSamp,1)*sum(M~=0);
%         Mor     = sum(M~=0)'+sum(M~=0) - Mand;
        overlap_narrow = overlap_narrow + Mand./Mor/experiment(useExpt).nreps;
    end

    subplot(2,2,1+2*(useExpt-1));
    plot(Cpn,Ckc,'.','color',[.75 .75 .75]);
    box off;hold on;axis equal
    plot(Cpn_broad,Ckc_broad,'r.');
    plot(Cpn_narrow,Ckc_narrow,'b.');
    title('corr')
    xlim([-.5 1]);
    ylim([-.25 1]);
    set(gca,'ytick',-.25:.25:1);
    set(gca,'xtick',-.5:.25:1);

    subplot(2,2,2+2*(useExpt-1));
    plot(Cpn,overlap,'.','color',[.75 .75 .75]);
    box off;hold on;axis equal
    plot(Cpn_broad,overlap_broad,'r.');
    plot(Cpn_narrow,overlap_narrow,'b.');
    title('IoU')
    xlim([-.5 1]);
    ylim([0 1]);
    set(gca,'ytick',0.:.25:1);
    set(gca,'xtick',-.5:.25:1);
    end