%% Mixing probaility table
clear; clc;
load C:\Users\jbland\Documents\LearningSims\MixSimplexGeneral_mixprobs
load C:\Users\jbland\Documents\LearningSims\MixBeta_mixprobs
load C:\Users\jbland\Documents\LearningSims\MixBetaLambdaConstant_mixprobs
load C:\Users\jbland\Documents\LearningSims\MixBetaRestricted_mixprobs
load C:\Users\jbland\Documents\LearningSims\MixBeta_Both_mixprobs


%%

% --- full table
%ColList = {'MixBeta','MixSimplexGeneral','MixBetaExpt2NoSymmetric','MixBetaExpt2short','MixBetaExpt2'};
%tablename = '';
%MP{1} = mix_MixBeta; MP{2} = mix_MixSimplexGeneral; MP{3} = mix_MixBetaExpt2NoSymmetric; MP{4} = mix_MixBetaExpt2short;
%MP{5} = mix_MixBetaExpt2;

% --- Just Expt 1 Beta
%ColList = {'MixBeta'};
%tablename = '_Expt1Beta';
%MP{1} = mix_MixBeta;

% --- Just Beta - Both experiments
ColList = {'MixBeta_Both'};
tablename = '_MixBeta_Both';
MP{1} = mix_MixBeta_Both;

% --- Just Expt 1 Beta LambdaConstant
%ColList = {'MixBetaLambdaConstant'};
%tablename = '_Expt1LambdaConstant';
%MP{1} = mix_MixBetaLambdaConstant;

% --- Just Expt 1 Beta Restricted
%ColList = {'MixBetaRestricted'};
%tablename = '_Expt1BetaRestricted';
%MP{1} = mix_MixBetaLambdaConstant;



% --- Just Expt 1 Simplex
%ColList = {'MixSimplexGeneral'};
%tablename = '_Expt1Simplex';
%MP{1} = mix_MixSimplexGeneral;

% --- Expt 2 Asym
%ColList = {'MixBetaExpt2NoSymmetric'};
%tablename = '_Expt2Asym';
%MP{1} = mix_MixBetaExpt2NoSymmetric;

% --- Expt 2 Short
%ColList = {'MixBetaExpt2short'};
%tablename = '_Expt2Short';
%MP{1} = mix_MixBetaExpt2short;

% --- Expt 2 Sym
%ColList = {'MixBetaExpt2'};
%tablename = '_Expt2Sym';
%MP{1} = mix_MixBetaExpt2;

%%


TypeList = {'AB CB','AB CM','AM CB','AM CM'};
TLx = ['AB CB';'AB CM';'AM CB';'AM CM']
TypeMarg = [2 4; 3 4];
TypeMargList = {'CM','AM'};




%%



fid = fopen(['figures/MixProbCompare' tablename '.tex'],'w');
fidM = fopen(['figures/MixProb_MarginalOnly' tablename '.tex'],'w');
fidlist = {fid,fidM};

for ff = 1:numel(fidlist)
fprintf(fidlist{ff},'\\begin{tabular}{l ccccc}\\hline\\hline\n');
if numel(ColList)>2
fprintf(fidlist{ff},'& \\multicolumn{2}{c}{Exp 1} & \\multicolumn{3}{c}{Exp 2}\\\\   \n');
fprintf(fidlist{ff},'& Asym & Simplex & Asym & Short & Sym \\\\ \\hline \n');
end
end
% mixing probabilities
fprintf(fid,'\\multicolumn{3}{l}{{\\sc Mixing probabilities -- joint}} \\\\ \n');
for tt = 1:numel(TypeList)
    fprintf(fid,'\\quad %s ',TypeList{tt});
       for ee = 1:numel(ColList)
           x = MP{ee}(:,tt);
           xm(ee) = mean(x);
           xs(ee) = std(x);
       end
    fprintf(fid,' & %1.3f',xm);
    if numel(ColList)>1
        fprintf(fid,'\\\\ \n');
    end
    fprintf(fid,' & (%1.3f)',xs);
    fprintf(fid,'\\\\ \n');
end

fprintf(fid,'\\hline');

for ff = 1:numel(fidlist)
fprintf(fidlist{ff},'\\multicolumn{3}{l}{{\\sc Mixing probabilities -- marginal}} \\\\ \n');
for tt = 1:numel(TypeMargList)
    fprintf(fidlist{ff},'\\quad %s ',TypeMargList{tt});
     for ee = 1:numel(ColList)
           x = sum(MP{ee}(:,TypeMarg(tt,:)),2);
           xm(ee) = mean(x);
           xs(ee) = std(x);
     end
     fprintf(fidlist{ff},' & %1.3f',xm);
     if numel(ColList)>1
        fprintf(fidlist{ff},'\\\\ \n');
     end
    fprintf(fidlist{ff},' & (%1.3f)',xs);
    fprintf(fidlist{ff},'\\\\ \n');
end
 fprintf(fidlist{ff},'\\quad Pr(%s $>$ %s)  ',TypeMargList{1},TypeMargList{2});
 for ee =1:numel(ColList)
     x = sum(MP{ee}(:,TypeMarg(1,:)),2)-sum(MP{ee}(:,TypeMarg(2,:)),2);
     fprintf(fidlist{ff},'& %1.3f',mean(x>0));
     
 end
fprintf(fidlist{ff},'\\\\ \n');
end
fprintf(fid,'\\hline');


fprintf(fid,'\\multicolumn{3}{l}{{\\sc Prob modal type}} \\\\ \n');
for tt = 1:numel(TypeList)
    fprintf(fid,'\\quad %s ',TypeList{tt});
       for ee = 1:numel(ColList)
           x = max(MP{ee},[],2)==MP{ee}(:,tt);
           xm(ee) = mean(x);
       end
    fprintf(fid,' & %1.3f',xm);
    fprintf(fid,'\\\\ \n');
end


for ff = 1:numel(fidlist)
fprintf(fidlist{ff},'\\hline\\hline');
fprintf(fidlist{ff},'\\end{tabular}\n');

    fclose(fidlist{ff});
end

%% table with 4! different ways to order the types

P = perms(1:4);
XM = zeros(size(P,1),numel(ColList));

fid = fopen(['figures/MixProbOrdering' tablename '.tex'],'w');

fprintf(fid,'\\begin{tabular}{l cccccc}\\hline\\hline\n');
fprintf(fid,'& \\multicolumn{2}{c}{Exp 1} & \\multicolumn{3}{c}{Exp 2}\\\\   \n');
fprintf(fid,'& Asym & Simplex & Asym & Short & Sym \\\\ \\hline \n');

for pp = 1:size(P,1)
    for qq = 1:3
        fprintf(fid,'%s $\\geq$ ',TypeList{P(pp,qq)});
    end
    fprintf(fid,'%s',TypeList{P(pp,4)});
    for ee = 1:numel(ColList)
        x = (MP{ee}(:,P(pp,1))>= MP{ee}(:,P(pp,2))) & (MP{ee}(:,P(pp,2))>= MP{ee}(:,P(pp,3))) & (MP{ee}(:,P(pp,3))>= MP{ee}(:,P(pp,4)));
         xm = mean(x);
         XM(pp,ee) = xm;
         %if xm>0.001
             fprintf(fid,' & %1.4f',xm);
         %else
         %    fprintf(fid,' & -');
         %end
    end
   
    fprintf(fid,'\\\\ \n');
end


fprintf(fid,'\\hline\\hline');
fprintf(fid,'\\end{tabular}\n');
fclose(fid);

%%

TableLength = 7; % 
[temp,II] = sort(XM(:,1));

fid = fopen(['figures/MixProbOrderingShort' tablename '.tex'],'w');


fprintf(fid,'\\begin{tabular}{l cccccc}\\hline\\hline\n');
if numel(MP)>1
    fprintf(fid,'& \\multicolumn{2}{c}{Exp 1} & \\multicolumn{3}{c}{Exp 2}\\\\   \n');
    fprintf(fid,'& Asym & Simplex & Asym & Short & Sym \\\\ \\hline \n');
else
    fprintf(fid,'Type ranking & Posterior probability \\\\ \\hline \n');
end

for tt = 1:TableLength
     pp = II(size(P,1)+1-tt);
     for qq = 1:3
        fprintf(fid,'%s $\\geq$ ',TypeList{P(pp,qq)});
     end
     fprintf(fid,'%s',TypeList{P(pp,4)});
     fprintf(fid,' & %1.4f',XM(pp,:));
     fprintf(fid,'\\\\ \n');
end

fprintf(fid,'\\hline\\hline');
fprintf(fid,'\\end{tabular}\n');
fclose(fid);

