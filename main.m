clear all
set(gcf,'visible','off') 

transition = [0.8, 0.2, 0.0; 0.1, 0.7, 0.2; 0.1, 0.0, 0.9];
emission = [0.1, 0.9; 0.5, 0.5; 0.9, 0.1];
pi = [0.7, 0.1, 0.2];

trans_guess = [0.799, 0.2, 0.001; 0.1, 0.7, 0.2; 0.1, 0.001, 0.899];
emission_guess = [0.4, 0.6; 0.5, 0.5; 0.6, 0.4];
pi_guess = [0.3, 0.3, 0.4];

N = size(transition,1);
K = size(emission,2);

given_runs = 1000;
random_runs = 10;
trans_all_given_BW = zeros(N, N, given_runs);
emis_all_given_BW = zeros(N, K, given_runs);
trans_all_random_BW = zeros(N, N, random_runs);
emis_all_random_BW = zeros(N, K, random_runs);
trans_all_given_GS = zeros(N, N, given_runs);
emis_all_given_GS = zeros(N, K, given_runs);
trans_all_random_GS = zeros(N, N, random_runs);
emis_all_random_GS = zeros(N, K, random_runs);
wrong_guess_given = zeros(3, given_runs); %first line viterbi, second line gibbs MAP, third line viterbi differs gibbs
wrong_guess_random = zeros(3, random_runs); %first line viterbi, second line gibbs MAP, third line viterbi differs gibbs

hmm_random = HMM(N, K);
hmm_given = HMM(N, K);
hmm = HMM(N, K);

length_vector = [100, 1000, 5000, 10000];
% length_vector = [5000];

for i = 1:size(length_vector, 2)
    length = length_vector(i);
    
    close all;       
    hold on;
    for run = 1:random_runs   
        disp(['random runs: ', num2str(run)]);
        %generate observations
        hmm = hmm.set(transition, emission, pi);
        [states, obser] = hmm.GenerateObservation(length); 
        
        hmm_random = HMM(N, K);
        %BW
        hmm_posterior = hmm_random;
        iter = 2000;        
        [hmm_posterior.transition, hmm_posterior.emission, hmm_posterior.pi] = hmm_posterior.BaumWelch(iter, length, obser);
        [guess_states_viterbi, guess_prob] = hmm_posterior.Viterbi(obser);
        wrong_ratio_viterbi = sum(states ~= guess_states_viterbi)/length;
        wrong_guess_random(1, run) = wrong_ratio_viterbi;
        trans_all_random_BW(:,:,run) = hmm_posterior.transition;    
        emis_all_random_BW(:,:,run) = hmm_posterior.emission;

        %GS
        hmm_posterior = hmm_random;
        iter = 200; % based on the estimated parameter, probability changed 0.005 per run
        [hmm_posterior.transition, hmm_posterior.emission, hmm_posterior.pi, guess_states_gibbs] = hmm_posterior.GibbsSampling(iter, obser);
        wrong_ratio_gibbs = sum(states ~= guess_states_gibbs)/length;
        wrong_guess_random(2, run) = wrong_ratio_gibbs;
        trans_all_random_GS(:,:,run) = hmm_posterior.transition;    
        emis_all_random_GS(:,:,run) = hmm_posterior.emission;
        
        differ_ratio_viterbi_gibbs = sum(guess_states_gibbs ~= guess_states_viterbi)/length;
        wrong_guess_random(3, run) = differ_ratio_viterbi_gibbs;
    end
    hold off;
    filepath = strcat('data/', num2str(length), '_BW_random.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)
    
    filepath = strcat('data/', num2str(length), '_BW_random_trans.txt');
    save(filepath, 'trans_all_random_BW');
    filepath = strcat('data/', num2str(length), '_BW_random_emis.txt');    
    save(filepath, 'emis_all_random_BW');    
    filepath = strcat('data/', num2str(length), '_GS_random_trans.txt');
    save(filepath, 'trans_all_random_GS');
    filepath = strcat('data/', num2str(length), '_GS_random_emis.txt');
    save(filepath, 'emis_all_random_GS');
    
    filepath = strcat('data/', num2str(length), '_random_wrong_guess.txt');
    save(filepath, 'wrong_guess_random');
    
    hold on;
    for run = 1:given_runs
        disp(['given runs: ', num2str(run)]);
        %generate observations
        hmm = hmm.set(transition, emission, pi);
        [states, obser] = hmm.GenerateObservation(length); 
        
        hmm_given = hmm_given.set(trans_guess, emission_guess, pi_guess);
        %BW
        hmm_posterior = hmm_given;
        iter = 2000;        
        [hmm_posterior.transition, hmm_posterior.emission, hmm_posterior.pi] = hmm_posterior.BaumWelch(iter, length, obser);
        [guess_states_viterbi, guess_prob] = hmm_posterior.Viterbi(obser);
        wrong_ratio_viterbi = sum(states ~= guess_states_viterbi)/length;
        wrong_guess_given(1, run) = wrong_ratio_viterbi;
        trans_all_given_BW(:,:,run) = hmm_posterior.transition;    
        emis_all_given_BW(:,:,run) = hmm_posterior.emission;

        %GS
        if(run <= 10)
            hmm_posterior = hmm_given;
            iter = 200; % based on the estimated parameter, probability changed 0.005 per run
            [hmm_posterior.transition, hmm_posterior.emission, hmm_posterior.pi, guess_states_gibbs] = hmm_posterior.GibbsSampling(iter, obser);
            wrong_ratio_gibbs = sum(states ~= guess_states_gibbs)/length;
            wrong_guess_given(2, run) = wrong_ratio_gibbs;
            trans_all_given_GS(:,:,run) = hmm_posterior.transition;    
            emis_all_given_GS(:,:,run) = hmm_posterior.emission;

            differ_ratio_viterbi_gibbs = sum(guess_states_gibbs ~= guess_states_viterbi)/length;
            wrong_guess_given(3, run) = differ_ratio_viterbi_gibbs;             
        end      
    end
    hold off;
    filepath = strcat('data/', num2str(length), '_BW_given.png');
    print(gcf, '-dpng', filepath);
    delete(gcf)
    
    filepath = strcat('data/', num2str(length), '_BW_given_trans.txt');
    save(filepath, 'trans_all_given_BW');
    filepath = strcat('data/', num2str(length), '_BW_given_emis.txt');    
    save(filepath, 'emis_all_given_BW');    
    filepath = strcat('data/', num2str(length), '_GS_given_trans.txt');
    save(filepath, 'trans_all_given_GS');
    filepath = strcat('data/', num2str(length), '_GS_given_emis.txt');
    save(filepath, 'emis_all_given_GS');
    
    filepath = strcat('data/', num2str(length), '_given_wrong_guess.txt');
    save(filepath, 'wrong_guess_given');    
end

% 
% % for i = 1:N
% %     for j = 1:N
% %         vars_trans(i,j) = var(trans_all(i, j, :));
% %     end
% % end
% % disp(['Variation 0f Transition Matrix: ']);
% % disp(vars_trans);
% % 
% % for i = 1:N
% %     for j = 1:K
% %         vars_emis(i,j) = var(emis_all(i, j, :));
% %     end
% % end
% % disp(['Variation 0f Emission Matrix: ']);
% % disp(vars_emis);
% 
% 
% 
