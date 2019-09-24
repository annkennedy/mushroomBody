

pList       = 1:30;
numTrials   = 10;
doPlot      = 0;


if(doPlot)
    figure(100);clf;hold on;
    h=plot(0,0,'k.-');
    h3=plot(0,0,'b.','markersize',15);
    h2=plot(0,0,'r.','markersize',15);
%     set(gca,'yscale','log');
end

realonly    = 1;
apl         = 1;
homeo       = 2;
shuffle     = 1;
uniform     = 1;
sp          = 1;
static      = 1;
mixtures    = 1;
nreps       = 5;
ncells = 2000;
ORN = buildORNs({}, opts.realonly);

figure(104);clf;hold on;
distractCount=0; colors='kbgr'; distractodors=[];
ystore={}; genstore = {};
for nDistract = [0 10 20 50]
    clear generalizationStats;
    for rep = 1:nreps
        rKC = KCmeanStore{realonly,apl,homeo,shuffle,uniform,sp,static,mixtures,rep};
        nonsilent = setdiff(1:110,find(sum(rKC(:,1:110)~=0)<10));
        rKC     = rKC(:,nonsilent);
        N       = length(nonsilent);
        w       = ones(1,ncells);
        w_n     = 20;  % alpha
        w_a     = 1;   % beta
        w_n_d   = 1;   % gamma
        w_a_d   = 1;   % delta

        alpha = 1e-5;
        clear y;
        for i = 1:length(pList)
            p = pList(i);
            disp([num2str(nDistract) ' | ' num2str(rep) ' | ' num2str(p)]);

            for trial = 1:numTrials
                d               = zeros(1,N);
                useodors        = randperm(N,p);
                genodors        = setdiff(1:N,useodors);
                if(nDistract)
                    distractodors = genodors(randperm(length(genodors),nDistract));
                else
                    distractodors = [];
                end
                d(useodors)     = 10;
                
                % initialize weights
                w = initialize_weights(rKC,w_n,w_a,distractodors);
                wInit=w;
                
                
                iters = 0; deltaw = 1;change=inf;
                while iters<3000*30 %&& ((((change>1e-6)&&(mod(iters,nDistract))) || (~mod(iters,nDistract))))
                    
                    iters   = iters+1;
                    frDistract = 1/3;
                    if(mod(iters,1/frDistract)==0)
                        odor = useodors(mod(ceil(iters*frDistract),length(useodors))+1);
                    else
                        if(nDistract)
                            odor = distractodors(randperm(nDistract,1));
                        else
                            odor = useodors(mod(ceil(iters*frDistract),length(useodors))+1);
                        end
                    end
                    
                    rON     = w * rKC(:,odor);
                    STDP    = w_n   * rKC(:,odor) - w_a   * rON * rKC(:,odor);
                    R_STDP  = w_n_d * rKC(:,odor) - w_a_d * rON * rKC(:,odor);

                    deltaw  = (STDP + R_STDP*d(odor))*alpha;
                    
                    wOld = w;
                    w = w + deltaw';
%                     w = max(w,0);
                    change = sum(abs(w-wOld));
                    
                    if(doPlot)
                        set(h,'xdata',1:size(rKC,2),'ydata',w*rKC); drawnow;
                        set(h2,'xdata',useodors,'ydata',w*rKC(:,useodors)); drawnow;
                        set(h3,'xdata',distractodors,'ydata',w*rKC(:,distractodors)); drawnow;
                    end
                end

                rONpre = wInit*rKC;
                rON = w*rKC;
                y(trial,i) = getBoundStats(rON, rONpre, useodors, distractodors, ORN.odorclass(nonsilent));
                generalizationStats((rep-1)*numTrials+trial,i) = y(trial,i).false_pos_count;

            end
        end

    end
    
    distractCount=distractCount+1;
    ystore{distractCount}   = y;
    genstore{distractCount} = generalizationStats;
    
    figure(104);hold on;
    drawvar(pList,generalizationStats,colors(distractCount),1)
    ylim([0 1])
    box off
    xlabel('# Paired Patterns','fontsize',14)
    ylabel('Probability of Generalization','fontsize',14);
    drawnow;
end
%%
for d = 1:4
    figure(104);hold on;
    drawvar(pList,genstore{d},colors(d),1)
    ylim([0 1])
    box off
    xlabel('# Paired Patterns','fontsize',14)
    ylabel('Probability of Generalization','fontsize',14);
    drawnow;
end

%%
for d = 1:4
    figure(104);hold on;
    drawvar(pList,genstore{d},colors(d),1)
    ylim([0 1])
    box off
    xlabel('# Paired Patterns','fontsize',14)
    ylabel('Probability of Generalization','fontsize',14);
    drawnow;
end
    
%%
mindist=zeros(30,1);
meandist=zeros(30,1);
meddist=zeros(30,1);
fpc = zeros(30,1);
for i=1:30
    for trial = 1:numTrials
        mindist(i)  = mindist(i)  + y(trial,i).min_distance/numTrials;
        meandist(i) = meandist(i) + y(trial,i).mean_distance/numTrials;
        meddist(i)  = meddist(i)  + y(trial,i).med_distance/numTrials;
        fpc(i) = fpc(i) + y(trial,i).false_pos_count/numTrials;
    end
end
figure(2);plot(mindist);
hold on;
plot(meandist,'g');
plot(meddist,'r');
plot([1 40],[0 0],'k--')

%%
samecat = zeros(30,1);
diffcat = zeros(30,1);
for i=1:30
    for trial = 1:numTrials
        samecat(i) = samecat(i) + mean(y(trial,i).false_pos_cats(y(trial,i).trainodor_cats))/numTrials;
        diffcat(i) = diffcat(i) + mean(y(trial,i).false_pos_cats(setdiff(1:10,y(trial,i).trainodor_cats)))/numTrials;
    end
end
figure(3);clf;
plot(samecat,'g');
hold on;
plot(diffcat);

%%

figure(100);clf;hold on;
color='rgbcmkrgbcmk';
for i=1:N
    plot(i,rON(i),[color(ORN.odorclass(nonsilent(i))) '.']);
end
plot(useodors,rON(useodors),'ro','markersize',15)
ylim([0 6])
xlim([0 N+1])