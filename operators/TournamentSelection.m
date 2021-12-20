function index = TournamentSelection(K,N,Fitness)
% TournamentSelection(The smaller the fitness value, the greater the probability of being selected.)
% Input: K round tournament, number of individuals to be selected, the fitness value.
% Output: selected individual index.
    Fitness = reshape(Fitness,1,[]);
    Fitness = Fitness + min(min(Fitness),0);
    [~,rank] = sort(1./Fitness);
    [~,rank] = sort(rank);
    Parents  = randi(length(1./Fitness),K,N);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*K);
end