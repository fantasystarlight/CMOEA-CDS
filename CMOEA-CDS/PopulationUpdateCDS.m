function Population = PopulationUpdateCDS(Population, N, k, constrained)
% The Constrained based search of CCMO-CDS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu
    %% Fitness calculation
    PopObj = Population.objs;
    [P, ~] = size(PopObj);
    if constrained
        PopCon = Population.cons;
        CV = sum(max(0,PopCon),2);
        Dominated = zeros(P);   
        CV(CV < 0) = 0;
    else
        CV = zeros(P,1);
    end
    crowdingdistance = pdist2(Population.objs, Population.objs);
    maxdistance = max(crowdingdistance,[],"all");
    for i = 1: P-1
        for j = i+1: P
            Domination = any(PopObj(i,:)>PopObj(j,:))-any(PopObj(i,:)<PopObj(j,:));
            if Domination == 1 && CV(i) >= CV(j)
                Dominated(i,j) = P;
            elseif Domination == -1 && CV(i) <= CV(j)
                Dominated(j,i) = P;
            elseif CV(i) > CV(j)
                Dominated(i,j) = 1 - (sqrt(sum((PopObj(i,:) - PopObj(j,:)).^2)) / maxdistance)^(1/k);
            elseif CV(i) < CV(j)
                Dominated(j,i) = 1 - (sqrt(sum((PopObj(i,:) - PopObj(j,:)).^2)) / maxdistance)^(1/k);
            else
                Dominated(i,j) = 1 - (sqrt(sum((PopObj(i,:) - PopObj(j,:)).^2)) / maxdistance)^(1/k);
                Dominated(j,i) = 1 - (sqrt(sum((PopObj(i,:) - PopObj(j,:)).^2)) / maxdistance)^(1/k);
            end
        end
    end
    PenaltyValue = sum(Dominated, 2);
    [~, Rank] = sort(PenaltyValue);
    Next(Rank(1:N)) = true;    
    Population = Population(Next);    
end