function [nWins, sumWins] = computeNWins(results)
    nMethods = size(results,2);
    nWins = nan(nMethods,nMethods);
    for method1=1:nMethods
        for method2=1:nMethods
            nWins(method1,method2) = sum(results(:,method2)>results(:,method1));
        end
    end
    sumWins = sum(nWins,1);
end

