clc
clear
dataset = {'Buettner','Darmanis','Engel','Pollen','Breton','Islet'};
for num = 1:6
    in_X=full(double(in_X));
    E = xishu(in_X);
    label=true_labs;
    [in_X,label]=paixu(in_X,label);
    n_space = length(unique(label));
    [X,] = FilterGenesZero(in_X);
    [n,m]=size(X);
    X = normalize_W(X', 2);
    Cdiss=cell(11,1);
    for d=5:5
        fprintf('Start Num %1.0f Datasets\n',num);
        for i=2:30
            options.maxIter = 50;
            [W,~] = NNDSVD(X,i,1);
            H = W'*X;
            [U,V,err] = NMF(X,options, i,W,H);
            Diss = corr(V, V, 'type' , 'Pearson');
            Diss = abs(Diss);
            Diss(Diss<(0.1*d))=0;
            Cdiss{i}= Diss;
            Cdiss1{i}= V;
            C1 = find(Diss==0);
            C2(i) = length(C1)/(n*n);
            C4(num,i) = C2(i);
            similarity2 = Diss+Diss';
            grps2 = SpectralClustering(similarity2, n_space);
            if i>=2
                C3(i) = (C2(i)-C2(i-1));
                Q(num,i) = C4(num,i)- C4(num,i-1);
                Q2(num,i) = Q(num,i)- Q(num,i-1);
                Q3(num,i) = abs(abs(Q2(num,i))-abs(Q2(num,i-1)));
                if Q3(num,i)<=E && i>2
                    ARI(num,d) = Cal_ARI(label, grps2);
                    NMI(num,d) = Cal_NMI(label, grps2);
                    K(num,d) = i;
                    fprintf('数据集ARI： %3.4f,数据集NMI： %3.4f, 阈值： %3.4f, 阈值2： %3.4f\n',ARI(num,d),NMI(num,d),Q3(num,i) ,E);
                    break
                end
            end
        end
    end
end