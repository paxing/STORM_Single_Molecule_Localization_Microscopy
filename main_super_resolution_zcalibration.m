clear all;
clc;
filename = 'Z_calibration.tif';

image_info=imfinfo(filename);
number_of_frames =  numel(image_info);

%%
figure()
for frame_number=1:number_of_frames
image_data=imread(filename, 'index', frame_number);
image_data = image_data - mean(image_data(:));
data(:,:,frame_number)=image_data;
imagesc(image_data);
axis image
title(frame_number)
colormap hot
pause(0.5)
end

%%
radius_of_the_ROI = 3;
size_of_the_ROI = radius_of_the_ROI*2+1;
threshold_value=1000;

localized_all = cell(number_of_frames,1);
figure('Renderer', 'painters', 'Position', [10 10 900 900])
for frame_number=1:number_of_frames
    disp(frame_number)
    tic
    
image_data=imread(filename, 'index', frame_number);

image_data = image_data - mean(image_data(:));
image_data_smooth = imgaussfilt(image_data,1);

image_data_dilate= imdilate(image_data_smooth,ones(size_of_the_ROI));
local_maxima=image_data_smooth==image_data_dilate;
above_threshold = image_data_smooth>threshold_value;

[row_peak, col_peak]=find(above_threshold.*local_maxima);
    ROI_size=20;

% if frame_number <30
%     ROI_size=20;
% elseif frame_number>=30 && frame_number<=120
%     ROI_size=10;
% 
% elseif frame_number>120
%     ROI_size = 20;
% end
pixel_centers = (1:1:ROI_size)-(ROI_size+1)/2;

number_of_peaks=numel(row_peak);
localized_frame = zeros(number_of_peaks,6);
for peak_number = 1: number_of_peaks
    %localization code
    
    ROI_index_x = pixel_centers+0.5+col_peak(peak_number);
    ROI_index_y = pixel_centers+0.5+ row_peak(peak_number);
    
    
    if ROI_index_x(1) <1
        ROI_index_x= ROI_index_x - ROI_index_x(1)+1;
        
    elseif ROI_index_x(end)>256
        ROI_index_x= ROI_index_x - (ROI_index_x(end)-256);
    end

     if ROI_index_y(1)<1
            ROI_index_y= ROI_index_y - ROI_index_y(1)+1;
            
     elseif ROI_index_y(end)>256
             ROI_index_y= ROI_index_y - (ROI_index_y(end)-256);
     end
    pixel_centers_x=pixel_centers+ROI_index_x((ROI_size)/2);
    pixel_centers_y=pixel_centers+ROI_index_y((ROI_size)/2);

    
    ROI  = double(image_data_smooth(ROI_index_y,ROI_index_x));
    
    starting_value_x = col_peak(peak_number);
    starting_value_y = row_peak(peak_number);
    stating_values_sigma_x=4;
    stating_values_sigma_y=4;
    starting_values_amplitude = max(ROI(:));
    starting_values_background = min(ROI(:));
    
    starting_values = [starting_value_x starting_value_y stating_values_sigma_x ...
        stating_values_sigma_y starting_values_amplitude starting_values_background];
    
    
    minimum_value_x =starting_value_x-stating_values_sigma_x;
    minimum_value_y =starting_value_y-stating_values_sigma_y;
    minimum_values_sigma_x = 0.25*stating_values_sigma_x;
    minimum_values_sigma_y = 0.25*stating_values_sigma_y;
    minimum_values_amplitude =0.25*starting_values_amplitude;
    minimum_values_background = 0.25*starting_values_background;
    
    minimum_values = [minimum_value_x minimum_value_y minimum_values_sigma_x ...
        minimum_values_sigma_y minimum_values_amplitude minimum_values_background];
    
    
    maximum_value_x = starting_value_x+stating_values_sigma_x;
    maximum_value_y = starting_value_y+stating_values_sigma_y;
    maximum_values_sigma_x = 4*stating_values_sigma_x;
    maximum_values_sigma_y = 4*stating_values_sigma_y;    
    maximum_values_amplitude = 4*starting_values_amplitude;
    maximum_values_background = min(ROI(:)) + 2*std(ROI(:));
    
    maximum_values = [maximum_value_x maximum_value_y maximum_values_sigma_x ...
        maximum_values_sigma_y maximum_values_amplitude maximum_values_background];

opts = optimset('Display','off');

[Localized] = lsqnonlin(@GAUSSIAN_FIT_2D, double(starting_values),...
    double(minimum_values), double(maximum_values),opts,...
    [pixel_centers_x;pixel_centers_y;ROI]);

localized_frame(peak_number,:)=Localized;

end
localized_all{frame_number}=localized_frame;

hold on
scatter(localized_frame(:,1),localized_frame(:,2),0.0001,'.w')
set(gca,'Color','k')

drawnow
toc
end

%%

figure()

for frame_number=1:number_of_frames
    hold on
    scatter(frame_number.*ones(size(localized_all{frame_number},1),1),localized_all{frame_number}(:,3),'.r')

end

xlabel('Framer number','Interpreter','Latex')
ylabel('$$\sigma_x$$ (pixels)','Interpreter','Latex')
xlim([1 number_of_frames])
ylim([0 6])

set(gca,'FontSize',12);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
print(gcf,'sigma_x_z_calibration','-dpng','-r300'); 

%%
figure()

for frame_number=1:number_of_frames
    hold on
    scatter(frame_number.*ones(size(localized_all{frame_number},1),1),localized_all{frame_number}(:,4),'.r')

end

xlabel('Framer number','Interpreter','Latex')
ylabel('$$\sigma_y$$ (pixels)','Interpreter','Latex')
xlim([1 number_of_frames])
ylim([0 5])

set(gca,'FontSize',12);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
print(gcf,'sigma_y_z_calibration','-dpng','-r300'); 

%%
figure()

for frame_number=1:number_of_frames
    hold on
    scatter(frame_number.*ones(size(localized_all{frame_number},1),1),...
        localized_all{frame_number}(:,3)-localized_all{frame_number}(:,4),'.r')

end

xlabel('Framer number','Interpreter','Latex')
ylabel('$$\Delta\sigma$$ (pixels)','Interpreter','Latex')
xlim([1 number_of_frames])
ylim([-3 3])

set(gca,'FontSize',12);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
print(gcf,'delta_sigma_z_calibration','-dpng','-r300'); 

%% linear fit

 all_frames =[];
 all_delta_sigma =[];
for frame_number=1:number_of_frames
     all_frames =[all_frames;frame_number*10/number_of_frames.*ones(size(localized_all{frame_number},1),1)];
 all_delta_sigma =[all_delta_sigma;localized_all{frame_number}(:,3)-localized_all{frame_number}(:,4)];
end
mdl = fitlm(all_frames,all_delta_sigma);

%%
predicted_values=predict(mdl,all_frames);
figure()

for frame_number=1:number_of_frames
    hold on
    scatter(frame_number*10/number_of_frames.*ones(size(localized_all{frame_number},1),1),...
        localized_all{frame_number}(:,3)-localized_all{frame_number}(:,4),'.r')

end
plot(all_frames,predicted_values, 'k','LineWidth',2)
text(0.5,2,'$$\Delta\sigma = 0.24223z-1.0272$$','Interpreter','Latex', 'FontSize',14)

xlabel('Z (nm)','Interpreter','Latex')
ylabel('$$\Delta\sigma$$ (pixels)','Interpreter','Latex')
xlim([0 10])
ylim([-3 3])

set(gca,'FontSize',12);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
print(gcf,'curve_z_calibration','-dpng','-r300'); 


%%



