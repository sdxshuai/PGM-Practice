clear all
set(gcf,'visible','off')
 

length_vector = [100, 1000, 5000, 10000];
total_runs = 1000;

N = 3;
K = 2;

trans_all = zeros(N, N, total_runs);
emis_all = zeros(N, K, total_runs);

for l = 1:size(length_vector, 2)
    length = length_vector(l);

    filepath = strcat('data/', num2str(length), '_BW_random_trans.txt');
    trans_all = load(filepath, '-mat');
    
    filepath = strcat('data/', num2str(length), '_BW_random_emis.txt');
    emis_all = load(filepath, '-mat');   
    
    trans_all = permute(trans_all.trans_all_random_BW, [2 1 3]);
    trans_all = reshape(trans_all, [9, 10]);
    emis_all = permute(emis_all.emis_all_random_BW, [2 1 3]);
    emis_all = reshape(emis_all, [6, 10]);

    hold on
    for k = 1:9
        plot(k, trans_all(k,:)', 'o');
    end
    hold off
    ylim([0 1])
    xlim([0 10])
    set(gca,'XTick',[0:1:10])
    filepath = strcat('data/', num2str(length), '_BW_random_trans_analysis.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)
    
    hold on
    for k = 1:6
        plot(k, emis_all(k,:)', 'o');
    end
    hold off
    ylim([0 1])
    xlim([0 7])
    set(gca,'XTick',[0:1:7])
    filepath = strcat('data/', num2str(length), '_BW_random_emis_analysis.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)   
end

