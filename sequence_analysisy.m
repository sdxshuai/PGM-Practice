clear all
set(gcf,'visible','off')
 

length_vector = [100, 1000, 5000, 10000];
total_runs = 1000;

N = 3;
K = 2;

vars_wrong_guess_given = zeros(1, 3);
vars_wrong_guess_random = zeros(1, 3);
mean_wrong_guess_given = zeros(1, 3);
mean_wrong_guess_random = zeros(1, 3);

wrong_guess_given = zeros(3, 1000); %first line viterbi, second line gibbs MAP, third line viterbi differs gibbs
wrong_guess_random = zeros(3, 10); %first line viterbi, second line gibbs MAP, third line viterbi differs gibbs

for l = 1:size(length_vector, 2)
    length = length_vector(l);
    
    filepath = strcat('data/', num2str(length), '_random_wrong_guess.txt');
    wrong_guess_random = load(filepath, '-mat');
    wrong_guess_random = wrong_guess_random.wrong_guess_random;
    filepath = strcat('data/', num2str(length), '_given_wrong_guess.txt');
    wrong_guess_given = load(filepath, '-mat');
    wrong_guess_given = wrong_guess_given.wrong_guess_given;
    
    for i = 1:3
        mean_wrong_guess_random(i) = mean(wrong_guess_random(i,:));
    end
    disp(['Mean 0f wrong_guess_random: ', num2str(length)]);
    disp(mean_wrong_guess_random);
    for i = 1:3
        vars_wrong_guess_random(i) = var(wrong_guess_random(i,:));
    end
    disp(['Variation 0f wrong_guess_random: ', num2str(length)]);
    disp(vars_wrong_guess_random);
end
for l = 1:size(length_vector, 2)
    length = length_vector(l);
    
    filepath = strcat('data/', num2str(length), '_random_wrong_guess.txt');
    wrong_guess_random = load(filepath, '-mat');
    wrong_guess_random = wrong_guess_random.wrong_guess_random;
    filepath = strcat('data/', num2str(length), '_given_wrong_guess.txt');
    wrong_guess_given = load(filepath, '-mat');
    wrong_guess_given = wrong_guess_given.wrong_guess_given;    
    
    for i = 1:3
        if(i == 1)
            mean_wrong_guess_given(i) = mean(wrong_guess_given(i,:));        
        else
            mean_wrong_guess_given(i) = mean(wrong_guess_given(i,1:10)); 
        end        
    end
    disp(['Mean 0f wrong_guess_given: ', num2str(length)]);
    disp(mean_wrong_guess_given);
    for i = 1:3
        if(i == 1)
            vars_wrong_guess_given(i) = var(wrong_guess_given(i,:));
        else
            vars_wrong_guess_given(i) = var(wrong_guess_given(i,1:10));
        end
    end
    disp(['Variation 0f wrong_guess_given: ', num2str(length)]);
    disp(vars_wrong_guess_given);
end