% % Define the original matrix A and create a range of perturbation magnitudes
% A = [1, 2; 3, 4]; % Original matrix
% num_samples = 100; % Number of samples
% num_perturbations = 50; % Number of perturbations
% perturbation_magnitudes = linspace(0, 0.01, num_perturbations); % Vary the perturbation magnitude
% 
% % Initialize arrays to store the results
% average_errors = zeros(1, num_perturbations);
% 
% for i = 1:num_perturbations
%     % Generate the perturbation matrix ΔA as a small random perturbation of A
%     DeltaA = perturbation_magnitudes(i) * randn(size(A));
% 
%     % Create a set of random domain vectors
%     domain_vectors = randn(size(A, 2), num_samples);
% 
%     % Compute the image of A and A + ΔA for each domain vector
%     image_A = A * domain_vectors;
%     image_A_plus_DeltaA = (A + DeltaA) * domain_vectors;
% 
%     % Compute the error in the image for each domain vector
%     error_in_image = image_A - image_A_plus_DeltaA;
% 
%     % Compute the norm (magnitude) of the error vectors
%     error_norms = vecnorm(error_in_image);
% 
%     % Average error over all domain vectors
%     average_errors(i) = mean(error_norms);
% end
% 
% % Create a plot to visualize the convergence of the error
% figure;
% plot(perturbation_magnitudes, average_errors, '-o');
% xlabel('Perturbation Magnitude');
% ylabel('Average Error in the Image');
% title('Convergence of Error in the Image as Perturbation Magnitude Decreases');
% 
% % Optionally, you can save the plot to a file using saveas or print commands
% % saveas(gcf, 'error_convergence_plot.png');
% % print('error_convergence_plot.pdf', '-dpdf');

% Define the original matrix A
A = [1, 1; 1, 1]; % Original matrix
num_samples = 100; % Number of samples
num_steps = 5000; % Number of time steps
perturbation_magnitude = 0.01; % Fixed perturbation magnitude

% Initialize arrays to store the results
time = 1:num_steps;
average_errors = zeros(1, num_steps);

for t = 1:num_steps
    % Generate a small random perturbation matrix ΔA at each time step
    DeltaA = perturbation_magnitude * randn(size(A));
    
    % Create a set of random domain vectors
    domain_vectors = randn(size(A, 2), num_samples);
    
    % Compute the image of A and A + ΔA for each domain vector
    image_A = A * domain_vectors;
    image_A_plus_DeltaA = (A + DeltaA) * domain_vectors;
    
    % Compute the error in the image for each domain vector
    error_in_image = image_A - image_A_plus_DeltaA;
    
    % Compute the norm (magnitude) of the error vectors
    error_norms = vecnorm(error_in_image);
    
    % Average error over all domain vectors
    average_errors(t) = mean(error_norms);
end

% Create a plot to visualize the convergence of the error over time
figure;
plot(time, average_errors, '-o');
xlabel('Time Step');
ylabel('Average Error in the Image');
title('Convergence of Error in the Image Over Time');

% Optionally, you can save the plot to a file using saveas or print commands
% saveas(gcf, 'error_convergence_plot.png');
% print('error_convergence_plot.pdf', '-dpdf');
