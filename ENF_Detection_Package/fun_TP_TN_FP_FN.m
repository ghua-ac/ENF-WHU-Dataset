function [TP,TN,FP,FN] = fun_TP_TN_FP_FN(result,ground_truth)
L = length(result);
TP = 0;
FP = 0;
TN = 0;
FN = 0;
for i = 1:L
    if result(i)==1 && ground_truth(i)==1
        TP = TP + 1;
    end
    if result(i)==1 && ground_truth(i)==0
        FP = FP + 1;
    end
    if result(i)==0 && ground_truth(i)==0
        TN = TN + 1;
    end
    if result(i)==0 && ground_truth(i)==1
        FN = FN + 1;
    end
end
end