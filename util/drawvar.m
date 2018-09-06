function h2=drawvar(X,Y,color,numDev,varargin)
if(isempty(varargin))
    h=area(X,[nanmean(Y,1)-nanstd(Y,1)*numDev; nanstd(Y,1)*numDev; nanstd(Y,1)*numDev]');
else
    h=area(X,[Y - varargin{1}*numDev; varargin{1}*numDev; varargin{1}*numDev]'); %lets you pass variance
end
    
    set(h,'FaceColor','none')
    set(h,'LineStyle','none')
    
    switch color
        case 'r'
            set(h([2 3]),'FaceColor',[1 .8 .8]);
        case 'b'
            set(h([2 3]),'FaceColor',[.8 .9 1]);
        case 'g'
            set(h([2 3]),'FaceColor',[.8 1 .8]);
        case 'c'
            set(h([2 3]),'FaceColor',[.8 1 1]);
        case 'm'
            set(h([2 3]),'FaceColor',[1 .8 .9]);
        case 'y'
            set(h([2 3]),'FaceColor',[1 1 .8]);
        case 'k'
            set(h([2 3]),'FaceColor',[.8 .8 .8]);
    end
    
    %now add the average overtop:
    hold on;
    h2=plot(X,nanmean(Y,1),['.-' color],'linewidth',1)
%     plot(X,nanmean(Y,1)/max(abs(nanmean(Y,1))),color,'linewidth',2)