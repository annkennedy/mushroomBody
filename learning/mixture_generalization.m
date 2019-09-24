
realonly = 1; useraw = 1; oldPNmodel = 0;
clc;


% mixlist = {'geranyl acetate','methyl salicylate','isopentyl acetate'};
mixlist = {ORN.odornames{public} ORN.odornames{private}};
mixfracs = .05:.05:.95;
initialize_parameters;
odors   = getconcinds(ORN,'all');
odorids = ORN.odornames(odors);

KC.inhtype  = 'divisive';
KC.inhflat  = 100; %for divisive inhibition
KC.sc       = 12.5;

KC = KCsim(KC,PN);

%%
for i = 1:length(mixlist)
    temp = strcmp(ORN.odornames,mixlist{i});
    if(sum(temp)>0)
        trainList{i} = find(temp);
        testing = find(1-cellfun(@isempty,strfind(ORN.odornames,mixlist{i})));
        mixes   = find(1-cellfun(@isempty,strfind(ORN.odornames,'mix')));
        testList{i}   = [trainList{i} intersect(testing,mixes)]; %toss out things that aren't the mixtures we made (concentration data, etc)
% testList{i} = 1:201;
    end
end





[generalizationStats astore] = olfPerceptron(trainList,testList,PN,KC);

%%
figure(1);clf;
hold on
colors='rgb';
for useodor=1:3
    half=size(astore{useodor},2)-1;
    if(useodor==1)
        ran1 = half/2:-1:2;
        ran2 = half:-1:half/2+2;
    elseif(useodor==2)
        ran1 = 2:half/2;
        ran2 = half:-1:half/2+2;
    else
        ran1 = 2:half/2;
        ran2 = half/2+2:half;
    end
    h(useodor)=plot(1:half/2,[1 mean(astore{useodor}(:,ran1))],colors(useodor),'linewidth',1);
    hold on;
    plot(1:half/2,[1 mean(astore{useodor}(:,ran2))],colors(useodor),'linewidth',1);
    xlim([1 half/2])
    set(gca,'xtick',[1 2:2:half/2]);
    set(gca,'xticklabel',{'','90:10','80:20','70:30','60:40','50:50','40:60','30:70','20:80','10:90'});
    set(gca,'fontSize',10)
    ylabel('Probability of generalization')
    xlabel('Mix (trained odor : distrator odor)')
end
legend(h,mixlist,'fontsize',10);

%%
clear genplots genplotstemp t_half;
trainpub_testpub   = [];
trainpub_testpriv  = [];
trainpriv_testpub  = [];
trainpriv_testpriv = [];

for trainodor = 1:20
    for j = 0.95:-0.05:0.05
        teststr             = ['mix ' num2str(j) ' ' ORN.odornames{trainList{trainodor}}];        
        hits                = find(1-cellfun(@isempty,strfind(ORN.odornames(testList{trainodor}),teststr)));
        
        redundant                   = cellfun(@length,strfind(ORN.odornames(testList{trainodor}(hits)),ORN.odornames{trainList{trainodor}}));
        hits(redundant>1)           = [];
        genplotstemp(round((1-j)/.05),:)    = mean(astore{trainodor}(:,hits));
        
%         publictest = strfind(ORN.odornames(testList{trainodor}(hits))
    end
    
    spacefind = strfind(ORN.odornames(testList{trainodor}(hits)),'0.95');
    for j=1:19
        testodors{j} = ORN.odornames{testList{trainodor}(hits(j))};
        testodors{j} = testodors{j}(spacefind{j}(end)+5:end);
    end
    
    genplots{trainodor} = zeros(19,length(trainList));
    for stim = 1:size(genplotstemp,2)
        temp = find(genplotstemp(:,stim)<.5,1);
        if(isempty(temp))
            temp = 19;
        end        
        
        stimind = find(1-cellfun(@isempty,strfind(ORN.odornames([trainList{:}]),testodors{stim})));
        
        t_half(trainodor,stimind) = temp;
        genplots{trainodor}(:,stimind) = genplotstemp(:,stim);
    end
    
    
    if(trainodor<=10)
        trainpub_testpub   = [trainpub_testpub  mean(genplots{trainodor}(:,setdiff(1:10,trainodor)),2)];
        trainpub_testpriv  = [trainpub_testpriv mean(genplots{trainodor}(:,11:20),2)];
    else
        trainpriv_testpub  = [trainpriv_testpub  mean(genplots{trainodor}(:,1:10),2)];
        trainpriv_testpriv = [trainpriv_testpriv mean(genplots{trainodor}(:,setdiff(11:20,trainodor)),2)];
    end
        
end

%%
figure(2);clf;
plot(mean(trainpub_testpub,2));
hold on;
plot(mean(trainpub_testpriv,2),'c');
plot(mean(trainpriv_testpub,2),'r');
plot(mean(trainpriv_testpriv,2),'m');
set(gca,'xtick',2:2:19);
set(gca,'xticklabel',{'90:10','80:20','70:30','60:40','50:50','40:60','30:70','20:80','10:90'});
legend('train pub distr pub','train pub distr priv','train priv distr pub','train priv distr priv');
        
%%
figure(1);clf;
subplot(2,1,1)
for stimnum=1:10
    plot(mean(genplots{stimnum}(:,setdiff(1:10,stimnum)),2),'g')
    hold on;
    plot(mean(genplots{stimnum}(:,11:20),2),'b')
end
set(gca,'xtick',2:2:19);
set(gca,'xticklabel',{'90:10','80:20','70:30','60:40','50:50','40:60','30:70','20:80','10:90'});
title('train on public odors');

subplot(2,1,2);
for stimnum=11:20
    plot(mean(genplots{stimnum}(:,1:10),2),'g')
    hold on;
    plot(mean(genplots{stimnum}(:,setdiff(11:20,stimnum)),2),'b')
end
set(gca,'xtick',2:2:19);
set(gca,'xticklabel',{'90:10','80:20','70:30','60:40','50:50','40:60','30:70','20:80','10:90'});
title('train on private odors');

        
        
        
        
        
        
        
        
        
        