
Nsamples = 10;
[broad, narrow] = find_public_and_private_odors(ORN,Nsamples);
numTrials = 1;
pList = 1;%:10;
runStats = true;
usemodel = true;
doBroad = false; doNarrow = false;

NKC = 2000;             %total number of KC's
sparseness = 0.1;       %use for making fake KCs
eta = 0.0025;            %learning rate
color = 'b';

expt = experiment(2);

count=0;
for nDistract = 0%[0 5 10 20 50]
    count=count+1;
    forgetStats = zeros(numTrials,length(pList));
    generalizationStats = ones(numTrials*10,length(pList),110);

    for i=1:110%length(pList)
        p = 1;%pList(i);
        fprintf('%d\n',i);
        
        for rep = 1:expt.nreps
            KC.mean_thresholded = expt.results{rep};
            nonsilent = setdiff(1:110,find(sum(KC.mean_thresholded(:,1:110))==0));
            if(~any(nonsilent==i))
                continue;
            end
            N = 110;%length(nonsilent);  %total patterns (don't count silent odors)
            [~,useBroad] = intersect(nonsilent,broad);
            [~,useNarrow] = intersect(nonsilent,narrow);

            for trial=1:numTrials

                t               = ones(1,110);%N);
                if(doBroad)
                    trainodors      = useBroad(randperm(Nsamples,p));
                elseif(doNarrow)
                    trainodors      = useNarrow(randperm(Nsamples,p));
                else
                    trainodors      = i;%randperm(N,p);
                end
                t(trainodors)   = -1;
                genodors        = setdiff(1:N,trainodors);
                distractodors   = genodors(randperm(N-p,nDistract));
                useodors        = [trainodors distractodors];
                genodors        = setdiff(1:N,useodors);

                if(usemodel)
                    phi         = KC.mean_thresholded(:,1:110);%(:,nonsilent);
                else %test generalization for various controls
%                     phi         = double(rand(NKC,N) <= sparseness);
                    phi = KC.mean_thresholded(:,nonsilent);
                    for odor=1:length(nonsilent)
                        phi(:,odor) = phi(randperm(NKC),odor);
                    end
                end

                k = ones(1,NKC); %discrimination vector
                a = sign(k*phi); %assigned category
                iters = 0;
                while (sum(a(useodors)~=t(useodors))>0) && (iters<50000)
                    iters   = iters+1;
                    iswrong = intersect(find(a~=t),useodors);
                    j       = iswrong(1);

                    kdel    = eta*t(j)*phi(:,j)';
                    k       = k+kdel;
                    a       = sign(k*phi); %assigned category
                end

                generalizationStats(trial+(rep-1)*numTrials,i,genodors) = a(genodors);
%                 generalizationStats(trial+(rep-1)*numTrials,i) = sum(a(genodors)~=t(genodors));
%                 generalizationStats(trial+(rep-1)*numTrials,i) = generalizationStats(trial+(rep-1)*numTrials,i)/(N-p); %probability
            end
        end
    end

    h(count)=drawvar(pList,generalizationStats,color,1);
    xlim([0 50])
    box off
    xlabel('number of learned odors')
    ylabel('probability of overgeneralization');    
    ylim([0 1]);
    xlim([0 10]);
    drawnow;
end
