function Offspring = GA(Parent,Task,j,select,Parameter)
%GA - Genetic operators for real, binary, and permutation based encodings.
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 1
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    %% Genetic operators for real encoding
    % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    % Polynomial mutation

    if select == 1
        Lower = zeros(2*N,D);
        Upper = ones(2*N,D);
    else
        Lower = repmat(Task.Lb(j,:),2*N,1);
        Upper = repmat(Task.Ub(j,:),2*N,1);
    end
    
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
% variable swap
    swap_indicator  = rand(N,D) < 0;
    temp1 = Offspring(1:floor(end/2),:);temp2 = Offspring(floor(end/2)+1:floor(end/2)*2,:);
    temp3 = temp2;
    temp3(swap_indicator) = temp1(swap_indicator);
    temp1(swap_indicator) = temp2(swap_indicator);
    Offspring(1:floor(end/2),:) = temp1;Offspring(floor(end/2)+1:floor(end/2)*2,:) = temp3;
    
end