function Population = EnvironmentalSelection(Population,Offspring,selection_process,N,M)
% EnvironmentalSelection
% Input:population chromes and objective values, offspring chromes and
% objective values, selection_process (option: elitist,roulette wheel),
% task number M.
% Output: population chromes and objective values.
%--------------------------------------------------------------------------
    nvec = [Population.rnvec(Population.flag == M,:);Offspring.rnvec(Offspring.flag == M,:)];
    fitness = [Population.factorial_costs(Population.flag == M,1);Offspring.factorial_costs(Offspring.flag == M,1)];
    if strcmp(selection_process,'elitist')
        [~,index]=sort(fitness);
    elseif strcmp(selection_process,'roulette wheel')
        index = RouletteWheelSelection(N,fitness);
    elseif strcmp(selection_process,'Tournament')
        index = TournamentSelection(2,N,fitness);
    end
    Population.rnvec(Population.flag == M,:) = nvec(index(1:N),:);
    Population.factorial_costs(Population.flag == M,1) = fitness(index(1:N),1);
end