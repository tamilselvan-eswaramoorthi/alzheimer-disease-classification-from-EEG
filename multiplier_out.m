function [full_collectData,bitSolution] = multiplier_out(divedentmulti,multiRate2,multiRate3)

if nargin > 2
    C1 = 0.2;
    C2 = 0.5;
end

accessData = [];
for k = 1:size(divedentmulti,1)
    mulData1 = bitshift(divedentmulti(k,:),2);
    mulData2 = bitshift(multiRate2(k,:),3);
    mulData3 = bitshift(multiRate3(k,:),4);
    for t = 1:length(mulData1)
        full_collectData(t,:) = mulData1(t)+mulData2(t)+mulData3(t);
    end
    accessData= [accessData full_collectData];
end
changeData = transpose(accessData);


for g = 1:length(changeData)
    for h = 1:size(changeData,2)
        if changeData(g,h) > zeros
            bitSolution(g,h) = ones;
        else
            bitSolution(g,h) = zeros;
        end
    end
end

