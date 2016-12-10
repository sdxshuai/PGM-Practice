classdef HMM
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        transition = [0.33, 0.33, 0.34; 0.33, 0.33, 0.34; 0.33, 0.33, 0.34];
        emission = [0.5, 0.5; 0.5, 0.5; 0.5, 0.5];
        pi = [1, 0, 0];
        N = 3;
        K = 2;
    end
    
    methods      
        
        function self = HMM(N, K)
            A = rand([N, N]);
            self.transition = A./repmat(sum(A, 2), [1, N]);
            A = rand([N, K]);
            self.emission = A./repmat(sum(A, 2), [1, K]);
            A = rand([1, N]);
            self.pi = A./repmat(sum(A, 2), [1, N]);
            self.N = N;
            self.K = K;
        end
        
        function self = set(self, transition, emission, pi)
            self.transition = transition;
            self.emission = emission;
            self.pi = pi;
            self.N = size(transition, 1);
            self.K = size(emission, 2);
        end
        
        function PrintHMM(self)
            disp('======================================================');
            disp('HMM transition matrix: ');
            disp(self.transition);
            disp('HMM emission matrix: ');
            disp(self.emission);
            disp('HMM state prob: ');
            disp(self.pi)
        end
        
        function [states, observations] = GenerateObservation(self, length)
            observations_rand = rand([1, length]);
            states_rand = rand([1, length]);
            observations = zeros(1, length);
            states = zeros(1, length);
            states(1) = 1;
            trans_accu = cumsum(self.transition, 2);
            obser_accu = cumsum(self.emission, 2);
            pi_accu = cumsum(self.pi, 2);
            states(1) = nnz(states_rand(1) > pi_accu) + 1;
            for i = 1:length   
                observations(i) = nnz(observations_rand(i) > obser_accu(states(i), :)) + 1;
                if(i < length)
                    states(i+1) = nnz(states_rand(i) > trans_accu(states(i),:)) + 1; 
                end               
            end
        end
        
        function [alpha, scales, prob] = Forward(self, length, observations)
            alpha = zeros(length, self.N);
            scales = zeros(1, length);
            
            alpha(1,:) = self.pi .* self.emission(:,observations(1))';
            scales(1) = sum(alpha(1,:));
            alpha(1,:) = alpha(1,:) ./ scales(1);
            
            for i = 2:length
                scales(i) = 0;
                alpha(i,:) = (alpha(i-1,:) * self.transition) .* self.emission(:,observations(i))';
                scales(i) = sum(alpha(i,:));
                alpha(i,:) = alpha(i,:) ./ scales(i);
            end
            
            prob = sum(log(scales));
        end
        
        function [beta, scales, prob] =  Backward(self, length, observations)
            beta = zeros(length, self.N);
            scales = zeros(1, length);
            
            beta(length,:) = 1;
            scales(length) = sum(beta(length,:));
            beta(length,:) = beta(length,:) ./ scales(length);
            
            for i = length-1:-1:1
                scales(i) = 0;
                beta(i,:) = (self.transition * (self.emission(:, observations(i+1)) .* (beta(i+1,:)')))';
                scales(i) = sum(beta(i,:));
                beta(i,:) = beta(i,:) ./ scales(i);
            end
            beta_zero = self.pi * (self.emission(:, observations(1)) .* (beta(1,:)'));
            prob = log(beta_zero) + sum(log(scales));            
        end
        
        function Xi = ComputeXi(self, length, observations, alpha, beta)
            Xi = zeros(self.N, self.N, length-1);
            for t = 1:length-1
                %{
                sum = (alpha(i,:) * (self.transition')) * ((beta(i+1,:)') .* self.emission(:, observations(i+1)));
                Xi(:,:, i) = repmat(alpha(i,:)', 1, self.N) .* self.transition .* repmat(self.emission(:, observations(i+1))' .* beta(i+1,:), self.N, 1);
                Xi(:,:, i) = Xi(:,:, i) ./ sum;
                %}
                sum = 0.0;
                for i = 1:self.N
                    for j = 1:self.N
                        Xi(i, j, t) = alpha(t, i) * self.transition(i, j) * self.emission(j, observations(t+1)) * beta(t+1, j); 
                        sum = sum + Xi(i, j, t);
                    end
                end
                Xi(:,:,t) = Xi(:,:,t) ./ sum;
            end
        end
        
        function Gamma = ComputeGamma(Xi)            
            Gamma = squeeze(sum(Xi, 2))';
        end
        
        function [trans, emis, pi] = BaumWelch(self, iter_num, length, observations)
            DELTA = 0.01;
            ratio = 0.0;
            delta = 0.0;
            delta_previous = 0.0;
            prob = 0.0;
            prob_previous = 0.0;
            MSP = 0.001;
            MOP = 0.001;
            ExitError = 10e-5;

            deltas = zeros(1, iter_num);
            for l = 1:iter_num
                [alpha, scales_alpha, prob] = self.Forward(length, observations);
                [beta, scales_beta, p_beta] = self.Backward(length, observations);
                if(prob - p_beta > 0.00000001)
                    break;
                end
                Xi = self.ComputeXi(length, observations, alpha, beta);
                % Gamma = self.ComputeGamma(Xi);
                Gamma = squeeze(sum(Xi, 2))';
                state_expect = sum(Gamma, 1);
                trans_expect = squeeze(sum(Xi, 3));
                trans_expect = MSP + (1 - self.N * MSP) .* trans_expect;
                
                self.pi = MSP + (1 - self.N * MSP) .* Gamma(1,:);
                self.transition = trans_expect ./ repmat(state_expect', 1, self.N);
                
                state_emission = zeros(self.N, self.K);
                for k = 1:self.K
                    obs_temp = (observations == k);
                    gamma_temp = Gamma .* repmat(obs_temp(1:length-1)', 1, self.N);
                    state_emission(:,k) = sum(gamma_temp, 1)';                    
                end
                state_emission = MOP + (1 - self.K * MOP) .* state_emission;
                self.emission =  state_emission ./ repmat(state_expect', 1, self.K);
                
                if(l == 1)
                    prob_previous = prob;
                    continue;
                end
                deltas(l) = prob - prob_previous;
                if(l == 2)
                    continue;
                end
                
                ratio = abs(deltas(l) / deltas(l - 1));
                prob_previous = prob;
                %display(['iter: ',num2str(l), '-----delta: ', num2str(deltas(l)), '-----ratio: ', num2str(ratio)]);
                if(abs(deltas(l)) < ExitError)
                    break;
                end
            end
            trans = self.transition;
            emis = self.emission;
            pi = self.pi;
            plot(log10(abs(deltas(1:end))));
        end
        
        function [states, prob] = Viterbi(self, observations)
            length = size(observations, 2);
            delta = zeros(length, self.N);
            phi = zeros(length, self.N);
            states = zeros(1, length);
            scales = zeros(1, length);     
            
            delta(1,:) = self.pi .* self.emission(:, observations(1))';
            scales(1) = sum(delta(1,:));
            delta(1,:) = delta(1,:) ./ scales(1);
            
            for t = 2:length
                for i = 1:self.N
                    delta_trans = delta(t-1,:) .* self.transition(:, i)';
                    [max_dt, idx] = max(delta_trans);
                    delta(t,i) = self.emission(i, observations(t)) * max_dt;
                    scales(t) = scales(t) + delta(t,i);
                    phi(t,i) = idx;
                end
                delta(t,:) = delta(t,:) ./ scales(t);
            end
            [prob, states(length)] = max(delta(length,:));
            prob = log(prob) + sum(log(scales));
            for t = length-1:-1:1
                states(t) = phi(t+1, states(t+1));
            end      
        end
        
        function [trans, emis, pi, states] = GibbsSampling(self, iter, observations)
            
            T = length(observations);
            %max_runs = ceil(400000 / T);
            max_runs = 200;
            states = zeros(1, T);
            
            for iteration = 1:iter
                y_guess = zeros(max_runs + 1, T);
                p = zeros(1, self.N); 
                trans = zeros(self.N, self.N);
                emis = zeros(self.N, self.K);
                pi = zeros(1, self.N);
                pi_accu = cumsum(self.pi, 2);        
                y_init_rand = rand([1, T]);
                for t = 1:T
                    y_guess(1, t) = nnz(y_init_rand(t) > pi_accu) + 1;
                end
                for n = 2:max_runs 
                    y_guess(n,:) = y_guess(n - 1,:);
                    for t = 1:T
                        if(t == 1)
                            for i = 1:self.N                       
                                p(i) = self.transition(i, y_guess(n, t+1)) * self.emission(i, observations(t));
                            end
                        elseif(t == T)
                            for i = 1:self.N                       
                                p(i) = self.emission(i, observations(t)) * self.transition(y_guess(n, t-1), i);
                            end
                        else
                            for i = 1:self.N                       
                                p(i) = self.transition(i, y_guess(n, t+1)) * self.emission(i, observations(t)) * self.transition(y_guess(n, t-1), i);
                            end
                        end

                        p = p./repmat(sum(p, 2), [1, self.N]);
                        p(isnan(p)) = 0;
                        
                        if(iteration == iter && n == max_runs)
                            [~, states(t)] = max(p);
                        end
                        p_accu = cumsum(p, 2);
                        y_rand = rand();
                        y_guess(n, t) = nnz(y_rand > p_accu) + 1;
                        if(y_guess(n, t) > 3)
                            y_guess(n, t) = 3;
                        end                        
                    end
                end

                start = 100;
                for n = start:max_runs
                    for t = 1:T
                        if (t ~= T)
                            trans(y_guess(n, t), y_guess(n, t+1)) = trans(y_guess(n, t), y_guess(n, t+1)) + 1;
                        end
                        emis(y_guess(n, t), observations(t)) = emis(y_guess(n, t), observations(t)) + 1;
                        pi(y_guess(n, t)) = pi(y_guess(n, t)) + 1;
                    end                
                end
                %smoothing
                trans(trans == 0) = 1;
                emis(emis == 0) = 1;
                pi(pi == 0) = 1;
                
                trans = trans ./ repmat(sum(trans, 2), [1, self.N]);
                emis = emis ./ repmat(sum(emis, 2), [1, self.K]);
                pi = pi ./ repmat(sum(pi, 2), [1, self.N]);
                
                self.transition = trans;
                self.emission = emis;
                self.pi = pi;
                disp(num2str(iteration));
                self.PrintHMM();
            end
            trans = self.transition;
            emis = self.emission;
            pi = self.pi;
        end
        
        function states = GibbsSamplingInfer(self, observations)
            T = length(observations);
            p = zeros(1, self.N);
            states = zeros(1, T);
            
            for t = 1:T
                if(t == 1)
                    for i = 1:self.N                       
                    	p(i) = self.transition(i, y_guess(n, t+1)) * self.emission(i, observations(t));
                    end
            	elseif(t == T)
                    for i = 1:self.N                       
                    	p(i) = self.emission(i, observations(t)) * self.transition(y_guess(n, t-1), i);
                    end
                else
                    for i = 1:self.N
                        p(i) = self.transition(i, y_guess(n, t+1)) * self.emission(i, observations(t)) * self.transition(y_guess(n, t-1), i);
                    end               
                end

            	p = p./repmat(sum(p, 2), [1, self.N]);
            	[~ ,states(t)] = max(p);
            end
        end
    end    
end

