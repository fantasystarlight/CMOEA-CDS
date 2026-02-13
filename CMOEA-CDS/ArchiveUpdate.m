function Population = ArchiveUpdate(Population, N)
% Archive update of CMOEA-CDS

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
    PopCon = Population.cons;
    CV = sum(max(0,PopCon),2);  
    CV(CV < 0) = 0;
    Next = false(1, length(Population));
    [~,Rank] = sort(CV);
    if sum(CV == 0) <= N
        Next(Rank(1:min(N,length(Population)))) = true;
        Population = Population(Next);
    else
        Next(Rank(1:sum(CV == 0))) = true;
        Population = Population(Next);
        [FrontNo, ~] = NDSort(Population.objs,N);
        Next = FrontNo == 1;
        if sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Population = Population(Next);
    end
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end