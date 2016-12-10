transition = [0.8, 0.2, 0.0; 0.1, 0.7, 0.2; 0.1, 0.0, 0.9];
emission = [0.1, 0.9; 0.5, 0.5; 0.9, 0.1];
pi = [0.7, 0.1, 0.2];

N = size(transition,1);
K = size(emission,2);

length_vector = [100, 1000, 5000, 10000];
total_runs = 1000;

trans_all = zeros(N, N, total_runs);
emis_all = zeros(N, K, total_runs);
pi_all = zeros(1, N, total_runs);

hmm = HMM(N, K);

for i = 1:size(length_vector, 2)
    length = length_vector(i);
    for run = 1:total_runs   
        disp(['runs--', num2str(length), '--', num2str(run)]);
        %generate observations
        hmm = hmm.set(transition, emission, pi);
        [states, obser] = hmm.GenerateObservation(length); 
        for t = 1:length
        	if (t ~= length)
                trans_all(states(t), states(t+1), run) = trans_all(states(t), states(t+1), run) + 1;
            end
            
        	emis_all(states(t), obser(t), run) = emis_all(states(t), obser(t), run) + 1;
        	pi_all(1, states(t), run) = pi_all(1, states(t), run) + 1;
        end
        trans_all(:,:,run) = trans_all(:,:,run) ./ repmat(sum(trans_all(:,:,run), 2), [1, N]);
        emis_all(:,:,run) = emis_all(:,:,run) ./ repmat(sum(emis_all(:,:,run), 2), [1, K]);
        pi_all(:,:, run) = pi_all(:,:, run) ./ repmat(sum(pi_all(:,:, run), 2), [1, N]);
    end
    filepath = strcat('data/', num2str(length), '_obsers_trans.txt');
    save(filepath, 'trans_all');
    filepath = strcat('data/', num2str(length), '_obsers_emis.txt');
    save(filepath, 'emis_all');
    filepath = strcat('data/', num2str(length), '_obsers_pi.txt');
    save(filepath, 'pi_all');
end