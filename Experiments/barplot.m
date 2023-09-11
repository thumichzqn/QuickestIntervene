clear all
load('data.mat')
Total_opt_mean = mean(Total_opt_lambda);
Total_Upper_mean = mean(Total_Upper_lambda);
Total_QCD_mean = mean(Total_QCD_lambda);
Total_opt_std = std(Total_opt_lambda);
Total_Upper_std = std(Total_Upper_lambda);
Total_QCD_std = std(Total_QCD_lambda);
figure;
hold on;
errorbar(viral_prob_set,Total_opt_mean,Total_opt_std)
errorbar(viral_prob_set,Total_Upper_mean,Total_Upper_std)
errorbar(viral_prob_set,Total_QCD_mean,Total_QCD_std)
set(gca,'Xscale','log')
legend('opt','Upper','QCD')