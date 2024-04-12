function ssim_value = calculate_ssim(I, I_prime, L, K1, K2)
    C1 = (K1 * L)^2;
    C2 = (K2 * L)^2;
    
    mu_x = mean2(I);
    mu_y = mean2(I_prime);
    sigma_x = std2(I);
    sigma_y = std2(I_prime);
    sigma_xy = sum(sum((I - mu_x) .* (I_prime - mu_y))) / (numel(I) - 1);
    
    numerator = (2 * mu_x * mu_y + C1) * (2 * sigma_xy + C2);
    denominator = (mu_x^2 + mu_y^2 + C1) * (sigma_x^2 + sigma_y^2 + C2);
    
    ssim_value = numerator / denominator;
end

