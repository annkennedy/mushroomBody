function [public, private] = find_public_and_private_odors(ORN,numodors)

if(~exist('numodors','var'))
    numodors=10;
end

temp = [(1:110)' sort(ORN.rates(ORN.HCList,1:110),'descend')'];
temp_srt = sortrows(temp,9);
public = temp_srt(end-numodors+1:end,1);

temp = sort(ORN.rates(ORN.HCList,1:110),'descend')';

temp = [(1:110)' temp(:,1)-temp(:,4) temp];
temp_srt = sortrows(temp,2);
private = temp_srt(end-numodors:end-1,1);