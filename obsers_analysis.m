clear all
set(gcf,'visible','off')
 

length_vector = [100, 1000, 5000, 10000];
total_runs = 1000;

N = 3;
K = 2;

trans_all = zeros(N, N, total_runs);
emis_all = zeros(N, K, total_runs);
pi_all = zeros(1, N, total_runs);

vars_trans = zeros(3, 3);
vars_emis = zeros(3, 2);
vars_pi = zeros(1, 3);

for l = 1:size(length_vector, 2)
    length = length_vector(l);

    filepath = strcat('data/', num2str(length), '_obsers_trans.txt');
    trans_all = load(filepath, '-mat');
    
    filepath = strcat('data/', num2str(length), '_obsers_emis.txt');
    emis_all = load(filepath, '-mat');
    
    filepath = strcat('data/', num2str(length), '_obsers_pi.txt');
    pi_all = load(filepath, '-mat');
    
    
    
%     for i = 1:N
%         for j = 1:N
%             vars_trans(i,j) = var(trans_all.trans_all(i, j, :)) .* 10^5;
%         end
%     end
%     disp(['Variation 0f Transition Matrix: ', num2str(length)]);
%     disp(vars_trans);
% 
%     
%     for i = 1:N
%         for j = 1:K
%             vars_emis(i,j) = var(emis_all.emis_all(i, j, :)) .* 10^5;
%         end
%     end
%     disp(['Variation 0f Emission Matrix: ', num2str(length)]);
%     disp(vars_emis);
% 
%     
%     for i = 1:N
%         vars_pi(i) = var(pi_all.pi_all(1, i, :)) .* 10^5;
%     end
%     disp(['Variation 0f Pi Matrix: ', num2str(length)]);
%     disp(vars_pi);
    

    trans_all = permute(trans_all.trans_all, [2 1 3]);
    trans_all = reshape(trans_all, [9, 1000]);
    emis_all = permute(emis_all.emis_all, [2 1 3]);
    emis_all = reshape(emis_all, [6, 1000]);
    pi_all = permute(pi_all.pi_all, [2 1 3]);
    pi_all = reshape(pi_all, [3, 1000]);

    boxplot(trans_all');
    ylim([0 1])
    filepath = strcat('data/', num2str(length), '_trans_analysis.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)
    
    boxplot(emis_all');
    ylim([0 1])
    filepath = strcat('data/', num2str(length), '_emis_analysis.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)    

    boxplot(pi_all');
    ylim([0 1])
    filepath = strcat('data/', num2str(length), '_pi_analysis.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)
end

