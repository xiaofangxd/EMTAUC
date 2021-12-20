function [objective,calls,EvBestFitness_evn,eval_cost] = CalObj(Task,rnvec,i,calls,EvBestFitness_evn,select,maxfes,eval_cost,costexp,N)
% Calculate the objective function value of the i-th task


    d = Task.Tdims(i);
    nvars = rnvec(:,1:d);
    NN = size(rnvec,1);
    if select == 1
        minrange = Task.Lb(i,1:d);
        maxrange = Task.Ub(i,1:d);
        y=repmat(maxrange-minrange,[NN,1]);
        x = y.*nvars + repmat(minrange,[NN,1]);
    else
        x = nvars;
    end
    objective=Task.fun(i).fnc(x);
    for j=1:NN 
%         if (mod(sum(calls),(maxfes-N*Task.M)/10) == (N*Task.M)) && (sum(calls)-N*Task.M)~=0
%             for k=1:Task.M
%                 TotalEvaluations_evn((sum(calls)-N*Task.M)/((maxfes-N*Task.M)/10),k) = EvBestFitness_evn(k);
%             end
%         end
%         if mod(sum(calls'.*((Task.p/Task.p(1)).^2)),costexp*N) == 0 && sum(calls'.*((Task.p/Task.p(1)).^2))/(costexp*N) ~= 0
%             eval_cost(sum(calls'.*((Task.p/Task.p(1)).^2))/(costexp*N)) = Task.fun(end).fnc(x(j,:));            
%         end
        if sum(calls'.*((Task.p/Task.p(1)).^2)) <= maxfes
            calls(i) = calls(i) + 1;
            if calls(i)==1
                EvBestFitness_evn(i) = objective(j);
            else
                EvBestFitness_evn(i) = min(EvBestFitness_evn(i),objective(j));
            end
        end
    end
end