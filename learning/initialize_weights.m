function w = initialize_weights(rKC,w_n,w_a,distractodors)
% initializes weights for LR_generalization function
% N       = size(rKC,2);
% rON     = ones(N,1)*(w_n/w_a)*3;
% w       = [lsqnonneg(rKC',rON)'];
% w       = [(pinv(rKC')*rON)'];
% w       = max(w + randn(size(w))/std(w)/50,0);

w=0.1*ones(1,2000);

if(~isempty(distractodors))
    alpha = 1e-3;

    for rep = 1:50000
        odor    = distractodors(randperm(length(distractodors),1));
        rON     = w*[rKC(:,odor)];
        STDP    = w_n * [rKC(:,odor)] - w_a * rON * [rKC(:,odor)];

        deltaw  = STDP*alpha;
        w = w + deltaw';
        w = max(w,0);
    end
end