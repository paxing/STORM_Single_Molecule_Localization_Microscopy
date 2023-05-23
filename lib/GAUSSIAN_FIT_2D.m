function [Delta] = GAUSSIAN_FIT_2D(guess,inputs)
pixel_centers_x = inputs(1,:);
pixel_centers_y = inputs(2,:);

pixel_values = inputs(3:end,:);
[pixel_centers_row pixel_centers_col] = ndgrid(pixel_centers_y,pixel_centers_x);

estimate_center = guess(1:2);
estimate_sigma_x = guess(3);
estimate_sigma_y = guess(4);
estimate_amplitude = guess(5);
estimate_background = guess(6);

Estimated_data = (1/(sqrt(estimate_sigma_x*estimate_sigma_y)*2*pi))*estimate_amplitude*exp(-(pixel_centers_col- ...
    estimate_center(1)).^2/(2*estimate_sigma_x.^2)).*exp(-(pixel_centers_row- ...
    estimate_center(2)).^2/(2*estimate_sigma_y.^2))+estimate_background;
Delta = pixel_values - Estimated_data;

end