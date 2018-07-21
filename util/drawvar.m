function drawvar(X,Y,color,nstd,varargin)
if(isempty(varargin))
    h=area(X,[mean(Y,1)-std(Y,1)*nstd; std(Y,1)*nstd; std(Y,1)*nstd]');
else
    h=area(X,[Y - varargin{1}*nstd; varargin{1}*nstd; varargin{1}*nstd]'); %lets you pass variance
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
    plot(X,mean(Y,1),['.-' color],'linewidth',1)
%     plot(X,mean(Y,1)/max(abs(mean(Y,1))),color,'linewidth',2)