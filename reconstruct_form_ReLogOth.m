cs_ratios = [10, 20, 30, 40, 50, 60];
% Image directory
image_dir = './dataset/Set11';
% Get list of images
image_files = dir(fullfile(image_dir, '*.tif'));

% Iterate over different compression ratios
for cs_ratio = cs_ratios
    matrix_dir = './sampling_matrix/PhReLogOth'; % Directory for storing measurement matrices

    % Initialize arrays for storing PSNR and SSIM values
    psnr_values = zeros(length(image_files), 1);
    ssim_values = zeros(length(image_files), 1);
    image_names = cell(length(image_files), 1);
    % Initialize variables for accumulating PSNR and SSIM
    total_psnr = 0;
    total_ssim = 0;

    sparsity_level = sparsity_levels(sparsity_index);
    tic; % Start timer
    % Iterate over all images
    for k = 1:length(image_files)
        image_name = fullfile(image_dir, image_files(k).name);
        I = imread(image_name);

        if size(I, 3) == 3
            I = rgb2gray(I);
        end

        % Divide the image into non-overlapping blocks
        NB = 128; % Block size
        [height, width] = size(I);
        blocks = mat2cell(I, repmat(NB, 1, floor(height/NB)), repmat(NB, 1, floor(width/NB)));
        
        % Load DFT measurement matrix
        filename = fullfile(matrix_dir, sprintf('Phase_ReLogistic_Oth_Sparse_DFT_Measurement_Matrix_%d_%d.mat', cs_ratio, NB));
        load(filename, 'Phi'); % Load measurement matrix
        
        % Sparse transformation and compressed sampling
        for i = 1:size(blocks, 1)
            for j = 1:size(blocks, 2)
                block = double(blocks{i, j});
                
                % DCT sparse transformation
                S = dct2(block);

                % Initialize the reconstructed block
                S_hat_block = zeros(size(S));

                % Apply Phi and OMP to recover each column of S
                for col = 1:size(S, 2)
                    Y = Phi * S(:, col); % Compressed sampling

                    % Extract real and imaginary parts
                    Y_real = real(Y);
                    Y_imag = imag(Y);

                    % Quantization of real and imaginary parts
                    R_min = min(Y_real(:));  % Minimum of real part
                    R_max = max(Y_real(:));  % Maximum of real part
                    I_min = min(Y_imag(:));  % Minimum of imaginary part
                    I_max = max(Y_imag(:));  % Maximum of imaginary part

                    % Actual quantization operation to get quantized data
                    Y_real_quantized = round(255 * (Y_real - R_min) / (R_max - R_min));
                    Y_imag_quantized = round(255 * (Y_imag - I_min) / (I_max - I_min));

                    % Convert to bit stream
                    bit_stream_real = de2bi(Y_real_quantized(:), 8, 'left-msb');
                    bit_stream_imag = de2bi(Y_imag_quantized(:), 8, 'left-msb');

                    % Combine bit streams
                    combined_bit_stream = [bit_stream_real; bit_stream_imag];% Extract bit streams of real and imaginary parts

                    bit_stream_real = combined_bit_stream(1:end/2, :);
                    bit_stream_imag = combined_bit_stream(end/2+1:end, :);

                    % Convert bit streams back to quantized values
                    Y_real_quantized = bi2de(bit_stream_real, 'left-msb');
                    Y_imag_quantized = bi2de(bit_stream_imag, 'left-msb');

                    % Dequantization operation to get original signal values
                    Y_real = (Y_real_quantized / 255 * (R_max - R_min)) + R_min;
                    Y_imag = (Y_imag_quantized / 255 * (I_max - I_min)) + I_min;

                    % Combine real and imaginary parts to recover complex signal
                    Y = Y_real + 1i * Y_imag;
                    % Set SL0 algorithm parameters
                    sigma_min = 0.01;  % Minimum sigma value
                    % Image reconstruction (using SL0)
                    [S_hat_col, iterCount] = SL0(Phi, Y, sigma_min);

                    S_hat_block(:, col) = S_hat_col; % Place the recovered column back in the block
                end

                blocks{i, j} = idct2(S_hat_block); % Inverse DCT transform

            end
        end

        % Reconstruct the image
        reconstructed_image = cell2mat(blocks);

        image_names[k] = image_files(k).name;

        % Magnitude of reconstructed image
        reconstructed_image_magnitude = abs(reconstructed_image);

        % I is the original image, reconstructed_image is the reconstructed image
        psnr_values(k) = calculate_psnr(double(I), double(reconstructed_image));
        ssim_values(k) = calculate_ssim(double(I), double(reconstructed_image), 255, 0.01, 0.03);
        total_psnr = total_psnr + psnr_values(k);
        total_ssim = total_ssim + ssim_values(k);

        output_filename = sprintf('SL0_%s_PSNR_%0.2f_SSIM_%0.4f_cs_ratio_%d_NB_%d.tif', image_files(k).name, psnr_values(k), ssim_values(k), cs_ratio, NB);
        imwrite(uint8(reconstructed_image_magnitude), fullfile('./result/gray/ReLogOthIRLS', output_filename));
    end

    elapsed_time = toc; % Stop timer and return elapsed time
    fprintf('Total elapsed time: %.2f seconds\n', elapsed_time);
    
    % Calculate average PSNR and SSIM
    average_psnr = total_psnr / length(image_files);
    average_ssim = total_ssim / length(image_files);

    % Output average PSNR and SSIM
    fprintf('cs_ratio: %d, Average PSNR/SSIM: %0.2f/%0.4f\n', cs_ratio, average_psnr, average_ssim);

    % Output results for all images
    for k = 1:length(image_files)
        fprintf('Image: %s, PSNR/SSIM: %0.2f/%0.4f\n', image_names[k], psnr_values(k), ssim_values(k));
    end
end
