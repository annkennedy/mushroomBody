function [public, private] = find_public_and_private_odors(ORN,numodors)

if(~exist('numodors','var'))
    numodors=10;
end

temp = [(1:110)' mean(ORN.rates(ORN.HCList,1:110),1)'];
temp_srt = sortrows(temp,2);
public = temp_srt(end-numodors+1:end,1);

temp = [(1:110)' (max(ORN.rates(ORN.HCList,1:110),[],1)./mean(ORN.rates(ORN.HCList,1:110),1))'];
temp_srt = sortrows(temp,2);
private = temp_srt(end-numodors:end-1,1);
temp = sort(ORN.rates(ORN.HCList,1:110),'descend')';
