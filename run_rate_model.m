realonly = 1;   % realonly == 1 only models glomeruli from which H&C recorded (n = 23),1
                % realonly == 0 makes fake data for the remainder of glomeruli (not recommended)

ORN = {}; %olfactory receptor neurons
ORN         = buildORN(ORN, realonly);

% play with ORN responses
% r2 = ORN.rates(ORN.HCList,1:110);
% r2 = reshape(r2(randperm(length(r2(:)))),23,110);
% r2 = bsxfun(@minus,r2,min(r2')');
% r2 = bsxfun(@times,r2,1./mean(r2')') * mean(r2(:));
% ORN.rates(ORN.HCList,1:110) = r2;

% new code to add mixtures of odors ---------------------------------------
% mixinds     = randperm(110,20);                  %select odors to make mixtures of
% ratios      = 0.1:0.1:0.9;                       %set mixing ratios
% mixlist     = ORN.odornames(mixinds);
% [ORN,pairs] = addmixtures(ORN, mixlist, ratios);
% -------------------------------------------------------------------------

nOdor = 110;%size(ORN.rates,2);
odorlist = 1:nOdor;

LN={};  %lateral neurons (inhbition between glomeruli)
LN          = buildLN(LN);
spRan = 0.1;%0.05:0.05:0.2;
% for use = 1
%     if(use)
%         KCmeanStore_inh = {};
%     else
%         KCmeanStore_none = {};
%     end
    for sp = 1:length(spRan)
        for rep = 1
            %set up PN and KC models
            [PN, KC] = initialize_cells(ORN,realonly);
            % uniform wPNKC
%             populated = ORN.HCList * PN.PN_per_glom;
%             populated = sort([populated-4 populated-3 populated-2 populated-1 populated]);
%             for i = 1:2000
%                 KC.wPNKC(populated,i) = KC.wPNKC(populated(randperm(length(populated))),i);
%             end            

            %first, measure the amount of spontaneous input to KCs
            [~,PN_t,~]  = getPNdynamics(ORN,PN,LN,1);
            PN_t = kron(PN_t,ones(PN.PN_per_glom,1));
            PN_spont    = mean(PN_t(:,(PN.tOn-PN.tmin)/PN.dt/2:(PN.tOn-PN.tmin)/PN.dt),2);
            KC.spont    = KC.wPNKC'*PN_spont;
            disp('Model initialized');

            %next, tune inhibition to set KC representation sparsity
            disp('Fitting APL inhibition (this step is slow)...');
            clear PN_sample PN_RW_sample;
            odorset = 2:3:110;     % code takes a while to run so I don't use all odors
%             odorset = randperm(110,60);
            sp_target = spRan(sp); % target sparsity of model KC population
            for i = odorset
                [~,PN_sample{i},~,PN_RW_sample{i}] = getPNdynamics(ORN,PN,LN,i);
            end
%             PN_t = PN_RW;
%             if(use)
%                 KC = setKCthreshold_bycell(PN,PN_sample, odorset, KC,sp_target*2);
%                 KC = fitSparseness(PN,PN_sample,odorset,KC,sp_target);
%             else
                KC.wInhKC = zeros(KC.ncells,1);
                KC.wKCInh = zeros(1,KC.ncells);
%                 KC = scalePNKCwts(PN, PN_sample, odorset, KC);
                KC = setKCthreshold(PN,PN_sample,odorset,KC,sp_target);
%             end
            disp('APL inhibition fit');


            %now, simulate responses to odors specified by odorlist
            odorspercell = zeros(KC.ncells,1);
            fr_active    = zeros(nOdor,1);
            KCmean_st    = zeros(KC.ncells,nOdor);

            for odorid = odorlist
                [~,PN_t,~,PN_RW]  = getPNdynamics(ORN,PN,LN,odorid);
%                 PN_t = PN_RW;
                PN_t = kron(PN_t,ones(PN.PN_per_glom,1));

                [KC_rasters,KC_t,~] = getKCratedynamics(PN_t,KC);

                KCmean = sum(KC_rasters,2);
                fr_active(odorid) = length(nonzeros(KCmean>0))/KC.ncells;

                odorspercell(KCmean>0) = odorspercell(KCmean>0) + 1;

                KCmean_st(:,odorid) = KCmean;
                if(mod(odorid,500)==0)
                    disp([num2str(odorid) '. ' ORN.odornames{odorid} ': ' num2str(fr_active(odorid)*100) '%']);
                end
            end
            KCmeanStore_scalePNwts_noinh{sp,rep} = KCmean_st;
%             if(use)
%                 KCmeanStore_inh{sp,rep} = KCmean_st;
%             else
%                 KCmeanStore_none{sp,rep} = KCmean_st;
%             end

            %save output for use with other functions (most of these I haven't included in this code)
        %     save(['KCsim_' num2str(rep) '.mat'],'PN','KC','fr_active','KCmean_st','ORN');
        end
    end
% end
% 
% s=svd(KCmean_st);
% dimensionality = sum(s)^2/sum(s.^2)

%% effects on similarity with inh

figure(1);clf;hold on;
i=2;
silentBars = [];
silentCellBars = [];

for use = {KCmeanStore_inh}
    vals = [];
    for rep = 1:3
        M = use{1}{i,rep}(:,1:110);
        vals = [vals pdist(M,'hamming')];
    end
    silentBars(end+1) = sum(sum(M~=0)==0);
    silentCellBars(end+1) = sum(sum(M'~=0)==0);
    h=plotCDFFit(vals);
end
set(h,'linewidth',2,'color','r');

% normalize ORN responses
M = KCmeanStore_normORNs{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% shuffle ORN responses
M = KCmeanStore_randORNs{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% static PN model
M = KCmeanStore_RW_PNs{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% uniform MB innervation by PNs
M = KCmeanStore_uniform_wPNKC{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% finer-tuned KC thresholds, no inh
M = KCmeanStore_fineThr{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% 10% of cells for each odor
for i=1:110
    vals = rand(1,2000);
    s = sort(vals);
    M(:,i) = vals>s(round(2000*.9));
end
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = sort(pdist(M,'hamming'));
h = plotCDFFit(vals);
set(h,'linewidth',2,'color','k');


legend('simple thr','norm ORNs','rand ORNs','static PN model','uniform wPNKC',...
        'fine-tuned KC thresholds',...
        'rand KCs','location','best');
xlabel('normalized hamming distance');
ylabel('fraction of odor pairs');
xlim([0 .5]);

figure(2);clf;
bar(silentBars);
set(gca,'xtick',1:7,'xticklabels',{'simple thr','norm ORNs','rand ORNs',...
    'static PN model','uniform wPNKC','fine-tuned KC thresholds'},'xticklabelrotation',30);
box off;

figure(4);clf;
bar(silentCellBars);
set(gca,'xtick',1:7,'xticklabels',{'simple thr','norm ORNs','rand ORNs',...
    'static PN model','uniform wPNKC','fine-tuned KC thresholds'},'xticklabelrotation',30);
box off;

%% effects on similarity, no inh

figure(1);clf;hold on;
i=2;
silentBars = [];
silentCellBars = [];

for use = {KCmeanStore_none}
    vals = [];
    for rep = 1:3
        M = use{1}{i,rep}(:,1:110);
        vals = [vals pdist(M,'hamming')];
    end
    silentBars(end+1) = sum(sum(M~=0)==0);
    silentCellBars(end+1) = sum(sum(M'~=0)==0);
    h=plotCDFFit(vals);
end
set(h,'linewidth',2,'color','r');

% normalize ORN responses
M = KCmeanStore_normORNs_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% shuffle ORN responses
M = KCmeanStore_randORNs_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% static PN model
M = KCmeanStore_RW_PNs_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% uniform MB innervation by PNs
M = KCmeanStore_uniform_wPNKC_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% finer-tuned KC thresholds, no inh
M = KCmeanStore_fineThr_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
plotCDFFit(vals);

% scale PN weights "presynaptically"
M = KCmeanStore_scalePNwts_noinh{1,1};
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = pdist(M,'hamming');
h=plotCDFFit(vals);
set(h,'linewidth',2);

% 10% of cells for each odor
for i=1:110
    vals = rand(1,2000);
    s = sort(vals);
    M(:,i) = vals>s(round(2000*.9));
end
silentBars(end+1) = sum(sum(M~=0)==0);
silentCellBars(end+1) = sum(sum(M'~=0)==0);
vals = sort(pdist(M,'hamming'));
h = plotCDFFit(vals);
set(h,'linewidth',2,'color','k');

legstr = {'simple thr','norm ORNs','rand ORNs','static PN model','uniform wPNKC',...
        'fine-tuned KC thresholds','presynaptic plasticity',...
        'rand KCs'};
legend(legstr{:},'location','best');
xlabel('normalized hamming distance');
ylabel('fraction of odor pairs');
xlim([0 .5]);

figure(2);clf;
bar(silentBars);
set(gca,'xtick',1:7,'xticklabels',legstr,'xticklabelrotation',30);
box off;

figure(3);clf;
bar(silentCellBars);
set(gca,'xtick',1:7,'xticklabels',legstr,'xticklabelrotation',30);
box off;

%% effects on similarity- looking at tuning KC thresholds

figure(1);clf;hold on;
i=2;

for use = {KCmeanStore_none}
    vals = [];
    for rep = 1:3
        M = use{1}{i,rep}(:,1:110);
        vals = [vals pdist(M,'hamming')];
    end
    h=plotCDFFit(vals);
    set(h,'linewidth',2,'color','r');
end

% finer-tuned KC thresholds, no inh
M = KCmeanStore_fineThr_noinh{1,1};
vals = pdist(M,'hamming');
plotCDFFit(vals);

% finer-tuned KC thresholds, no inh
for rep=1:10
M = KCmeanStore_fineThr_noinh_sample60{1,rep};
vals = pdist(M,'hamming');
plotCDFFit(vals);
drawnow;
end

% 10% of cells for each odor
for i=1:110
    vals = rand(1,2000);
    s = sort(vals);
    M(:,i) = vals>s(round(2000*.9));
end
vals = sort(pdist(M,'hamming'));
h=plotCDFFit(vals);
set(h,'linewidth',2,'color','k');

% legend('no inh','with inh',...
%         'fine-tuned KC thresholds, no inh','same, fewer thr examples',...
%         'rand KCs','location','best');
xlabel('normalized hamming distance');
ylabel('fraction of odor pairs');
xlim([0 .5]);
%% number of missed odors
figure(2);clf;hold on;
for rep=1:10
    M = KCmeanStore_fineThr_noinh_sample5{1,rep};
    silent(rep) = sum(sum(M~=0)==0);
end
bar(1,mean(silent));
errorbar(1,mean(silent),std(silent),'k');

for rep=1:10
    M = KCmeanStore_fineThr_noinh_sample20{1,rep};
    silent(rep) = sum(sum(M~=0)==0);
end
bar(2,mean(silent));
errorbar(2,mean(silent),std(silent),'k');
ylabel('number of silent odors');

for rep=1:10
    M = KCmeanStore_fineThr_noinh_sample60{1,rep};
    silent(rep) = sum(sum(M~=0)==0);
end
bar(3,mean(silent));
errorbar(3,mean(silent),std(silent),'k');
ylabel('number of silent odors');
set(gca,'xtick',1:3,'xticklabel',{'5','20','60'});

%%
[public, private] = find_public_and_private_odors(ORN);
M = KCmeanStore_fineThr_noinh{1,1};
figure(2);clf;
bar(sum(M(:,[public; private])~=0))

