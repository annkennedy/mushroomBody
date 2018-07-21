function [w gListScale] = buildConnects(KC,PN, cxntype)

gData3;

nK = KC.ncells;                                 %number of Kenyon cells
nP = PN.ncells*PN.PN_per_glom;                  %number of projection neurons
ngroups = 1; %a default value

if(cxntype==1)                      %only take input from the real glomeruli
    
    numGlomKeep = sum(ismember(gList,PN.HCList));
    gListKeep   = gList(ismember(gList,PN.HCList));
    
%do structured connections; decimal gives p(structure)
elseif((cxntype>1)&(cxntype<=2))

    nG = length(PN.HCList);
    ngroups = 5; %let's stick with this for now
    glomGroups = reshape([PN.HCList(randperm(nG,nG)) zeros(1,ceil(nG/ngroups)*ngroups - nG)],ngroups,[]);
    
    for i = 1:ngroups
        numGlomKeep(i) = sum(ismember(gList,glomGroups(i,:)));
        gListKeep(i,1:numGlomKeep(i)) = gList(ismember(gList,glomGroups(i,:)));
    end
    
    
%do altered KC claw numbers
elseif(cxntype>3)
    
    numGlomKeep = sum(ismember(gList,PN.HCList));
    gListKeep   = gList(ismember(gList,PN.HCList));
    
    nClaws = round((cxntype-3)*100) %encode claw number in the decimal, i suck
    clawsR = nClaws*ones(size(clawsR));
    
else                                %take input from everyone
    
    numGlomKeep = sum(ismember(gList,1:51));
    gListKeep   = gList(ismember(gList,1:51));
    
end
allCells = nonzeros(gListKeep')';

nBar = mean(clawsR);                            %clawsR = array of #claws (inputs) per KC

wP   = zeros(nK, nP);                           %PN -> KC weight matrix
for i=1:nK
    KCgrp = ceil(i/(nK/ngroups));
    isGlobalKC = (rand(1)>(cxntype-1));          %cxntype <= 1 ->all global, cxntype = 2 ->all local
    iC = randi(numKCsRecorded, 1, 1);
        nConnect = clawsR(iC);                      %pick a number of claws for model KC
    
    if(~isGlobalKC)                               %for current KC, pick input gloms
        iR = randi(numGlomKeep(KCgrp), nConnect, 1);
    else
        iR = randi(sum(numGlomKeep), nConnect, 1);
    end
    iG = randi(PN.PN_per_glom,nConnect,1);      %and then pick which PN it gets from said glom

    
    for j = 1:nConnect                          %account for multiple inputs from the same PN
        if(~isGlobalKC)
            ind     = (gListKeep(KCgrp,iR(j))-1) * PN.PN_per_glom + iG(j);
            count   = sum((gListKeep(KCgrp,iR)==gListKeep(KCgrp,iR(j)))'.*(iG==iG(j)));
        else
            ind     = (allCells(iR(j))-1) * PN.PN_per_glom + iG(j);
            count   = sum((allCells(iR)==allCells(iR(j)))'.*(iG==iG(j)));
        end
        wP(i, ind)  = nBar*count/nConnect;
    end
end

w = wP';