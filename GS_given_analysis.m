clear all
set(gcf,'visible','off')
 

length_vector = [100, 1000, 5000, 10000];
total_runs = 1000;

N = 3;
K = 2;

trans_all = zeros(N, N, total_runs);
emis_all = zeros(N, K, total_runs);

vars_trans = zeros(3, 3);
vars_emis = zeros(3, 2);

for l = 1:size(length_vector, 2)
    length = length_vector(l);

    filepath = strcat('data/', num2str(length), '_GS_given_trans.txt');
    trans_all = load(filepath, '-mat');
    
    filepath = strcat('data/', num2str(length), '_GS_given_emis.txt');
    emis_all = load(filepath, '-mat');
    
    for i = 1:N
        for j = 1:N
            vars_trans(i,j) = var(trans_all.trans_all_given_GS(i, j, 1:10)) .* 10^5;
        end
    end
    disp(['Variation 0f Transition Matrix: ', num2str(length)]);
    disp(vars_trans);

    
    for i = 1:N
        for j = 1:K
            vars_emis(i,j) = var(emis_all.emis_all_given_GS(i, j, 1:10)) .* 10^5;
        end
    end
    disp(['Variation 0f Emission Matrix: ', num2str(length)]);
    disp(vars_emis);


%     trans_all = permute(trans_all.trans_all_given_GS, [2 1 3]);
%     trans_all = reshape(trans_all, [9, 1000]);
%     emis_all = permute(emis_all.emis_all_given_GS, [2 1 3]);
%     emis_all = reshape(emis_all, [6, 1000]);
% 
%     hold on
%     for k = 1:9
%         plot(k, trans_all(k,1:10)', 'o');
%     end
%     hold off
%     ylim([0 1])
%     xlim([0 10])
%     set(gca,'XTick',[0:1:10])
%     filepath = strcat('data/', num2str(length), '_GS_given_trans_analysis.png');
%     print(gcf, '-dpng', filepath);
%     delete(gcf)
%     
%     hold on
%     for k = 1:6
%         plot(k, emis_all(k,1:10)', 'o');
%     end
%     hold off
%     ylim([0 1])
%     xlim([0 7])
%     set(gca,'XTick',[0:1:7])
%     filepath = strcat('data/', num2str(length), '_GS_given_emis_analysis.png');
%     print(gcf, '-dpng', filepath);
%     delete(gcf)   
end

