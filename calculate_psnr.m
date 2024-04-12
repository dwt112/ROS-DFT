% function psnr_value = calculate_psnr(I, I_prime, N)
%     mse = sum(sum((I - I_prime) .^ 2)) / (N ^ 2);
%     psnr_value = 10 * log10((256 ^ 2) / mse);
% end
function psnr_value = calculate_psnr(I, I_prime)
    [M, N] = size(I);
    mse = sum(sum((I - I_prime) .^ 2)) / (M * N);
    psnr_value = 10 * log10((255 ^ 2) / mse);
end
