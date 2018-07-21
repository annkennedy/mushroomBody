    realonly = 1;   % realonly == 1 only models glomeruli from which H&C recorded (n = 23),
                % realonly == 0 makes fake data for the remainder of glomeruli (not recommended)

ORN={}; %olfactory receptor neurons
ORN         = buildORN(ORN, realonly);
nOdor       = size(ORN.rates,2);

LN={};  %lateral neurons
LN          = buildLN(LN);

%set up PN, KC models
[PN KC] = initialize_cells(ORN,realonly);

%%

KC.wInhKC = ones(KC.ncells,1) * 20;
KC.wKCInh = (1/KC.ncells) * ones(1,KC.ncells) * 25;

KC.wInhKC = KC.wInhKC.*((base_rates' + 44)/50).^3;

temp = min(base_rates,6);
KC.wInhKC = KC.wInhKC.*((temp' + 4)/10).^2;

KCmean_st = run_rate_model(KC,PN,ORN,LN);
% base_rates = sum(KCmean_st'>0);
mean(sum(KCmean_st>0)/1000)


figure(1);clf;
subplot(2,1,1);
plot(sort(sum(KCmean_st'>0)))

subplot(2,1,2);
plot(sort(sum(KCmean_st>0)))
axis tight