clear all;
clc;
%filename = 'STORM_Data_LAST200.tif';
filename = 'sequences\Aquired_STORM';

%image_info=imfinfo(filename);

radius_of_the_ROI = 3;
size_of_the_ROI = radius_of_the_ROI*2+1;
threshold_value=2000;

number_of_frames = 25000;
localized_all = cell(number_of_frames,1);
figure('Renderer', 'painters', 'Position', [10 10 900 900])
for frame_number=1:number_of_frames
    disp(frame_number)
    tic
    
%image_data=imread(filename, 'index', frame_number);
image_data=imread([filename, num2str(10000+frame_number), '.tif']);

image_data = image_data - mean(image_data(:));
image_data_smooth = imgaussfilt(image_data,1);

image_data_dilate= imdilate(image_data_smooth,ones(size_of_the_ROI));
local_maxima=image_data_smooth==image_data_dilate;
above_threshold = image_data_smooth>threshold_value;

[row_peak, col_peak]=find(above_threshold.*local_maxima);

ROI_size=10; %for the gaussian fitting
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
    stating_values_sigma_x=2;
    stating_values_sigma_y=2;
    starting_values_amplitude = max(ROI(:));
    starting_values_background = min(ROI(:));
    
    starting_values = [starting_value_x starting_value_y stating_values_sigma_x ...
        stating_values_sigma_y starting_values_amplitude starting_values_background];
    
    
    minimum_value_x =starting_value_x-stating_values_sigma_x;
    minimum_value_y =starting_value_y-stating_values_sigma_y;
    minimum_values_sigma_x = 0.5*stating_values_sigma_x;
    minimum_values_sigma_y = 0.5*stating_values_sigma_y;
    minimum_values_amplitude =0.5*starting_values_amplitude;
    minimum_values_background = 0.25*starting_values_background;
    
    minimum_values = [minimum_value_x minimum_value_y minimum_values_sigma_x ...
        minimum_values_sigma_y minimum_values_amplitude minimum_values_background];
    
    
    maximum_value_x = starting_value_x+stating_values_sigma_x;
    maximum_value_y = starting_value_y+stating_values_sigma_y;
    maximum_values_sigma_x = 2*stating_values_sigma_x;
    maximum_values_sigma_y = 2*stating_values_sigma_y;    
    maximum_values_amplitude = 2*starting_values_amplitude;
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


%% all frame

all_x = [];
all_y = [];

for frame_number = 1:number_of_frames
    
    all_x = [all_x; localized_all{frame_number}(:,1)];
    all_y = [all_y; localized_all{frame_number}(:,2)];
    
end

%%
figure('Renderer', 'painters', 'Position', [10 10 900 900])
scatter(all_x, all_y, 0.0001,'.w')
xlim([0.5 255.5])
ylim([0.5 255.5])
xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'Color','k')

set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
ax.XColor = 'w'; % 
ax.YColor = 'w'; % 
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'InvertHardCopy', 'off'); 

set(gcf,'Color',[0 0 0]);
print(gcf,'localization_scatter_25000','-dpng','-r300'); 

%%
figure('Renderer', 'painters', 'Position', [10 10 900 900])
scatter(all_x, all_y,10, '.w')
xlim([60.5 100.5])
ylim([100.5 140.5])
xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'Color','k')

set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.XColor = 'w'; % 
ax.YColor = 'w'; % 
ax.LineWidth=1.5;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'InvertHardCopy', 'off'); 

set(gcf,'Color',[0 0 0]);
print(gcf,'localization_scatter_25000_magnification','-dpng','-r300'); 
%%

[N,Xedges,Yedges] = histcounts2(all_x,all_y,[1024 1024]);

figure('Renderer', 'painters', 'Position', [10 10 900 900])

imagesc([0.5 255.5],[0.5 255.5], N.')
caxis([0 5])
colormap gray
xlim([0.5 255.5])
ylim([0.5 255.5])
xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'Color','k')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
ax.XColor = 'w'; % 
ax.YColor = 'w'; % 
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'InvertHardCopy', 'off'); 
set(gcf,'Color',[0 0 0]);
print(gcf,'localization_histo_25000','-dpng','-r300'); 
%%

[N,Xedges,Yedges] = histcounts2(all_x,all_y,[1024*5 1024*5]);
N = conv2(N, ones(5)/5^2);
figure('Renderer', 'painters', 'Position', [10 10 900 900])
imagesc([0.5 255.5],[0.5 255.5], N.')
xlim([60.5 100.5])
ylim([100.5 140.5])
caxis([0.1 0.75])
colormap gray
xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'Color','k')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.XColor = 'w'; % 
ax.YColor = 'w'; % 
ax.LineWidth=1.5;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'InvertHardCopy', 'off'); 

set(gcf,'Color',[0 0 0]);
print(gcf,'localization_histo_25000_magnification','-dpng','-r300'); 

%% 3d
all_z=[];
for frame_number = 1:number_of_frames
    
    all_z = [all_z; 4.13*(localized_all{frame_number}(:,3)-localized_all{frame_number}(:,4))+4.24];
    
end

%%
zaxis=1:0.05:8.5;
n3d=histcnd(all_x, all_y,all_z, 0.5:0.25:255.5,0.5:0.25:255.5,zaxis);

n3d=convn(n3d,ones(5,5,3)./5^2*3);
figure('Position', [10 10 900 1600])
for i= 1:10:size(n3d,3)
    nexttile
    imagesc(n3d(:,:,i))
    title([num2str(zaxis(i)) ' nm'], 'Interpreter', 'Latex', 'FontSize',18)
    axis image
    colormap gray
    caxis([0 mean(n3d(:))+2*std(n3d(:))])
    set(gca,'YDir','reverse')
    %set(gca, 'visible', 'off')
    set(gca,'xtick',[])
set(gca,'xticklabel',[])

    set(gca,'ytick',[])
set(gca,'yticklabel',[])

end
print(gcf,'localization_histo_25000_3d','-dpng','-r300'); 

%%
h = slice(n3d, [], [], 1:5:size(n3d,3));
set(h, 'EdgeColor','none', 'FaceColor','interp')
colormap(flipud(gray))
caxis([0 mean(n3d(:))+std(n3d(:))])
alpha(.1)
%%
figure('Position', [10 10 900 900])

scatter3(all_x, all_y,all_z,0.001,'.w');

xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
zlabel('z(nm)','Interpreter','Latex')

set(gca,'Color','k')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
 
ax.LineWidth=1.5;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'InvertHardCopy', 'off');
%% single frame

figure('Renderer', 'painters', 'Position', [10 10 900 900])
hold on
imagesc([0.5 255.5],[0.5 255.5],image_data);
colormap gray
axis image
scatter(localized_frame(:,1),localized_frame(:,2),'.r')

xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
set(gca,'LooseInset',get(gca,'TightInset'));
print(gcf,'localization_frame200','-dpng','-r300'); 
%%
figure('Renderer', 'painters', 'Position', [10 10 900 900])
hold on
imagesc([0.5 255.5],[0.5 255.5],image_data);
colormap gray
axis image
scatter(localized_frame(:,1),localized_frame(:,2),'.r')
xlim([30.5 45.5])
ylim([115.5 130.5])
xlabel('Pixel coordinates','Interpreter','Latex')
ylabel('Pixel coordinates','Interpreter','Latex')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
set(gca,'LooseInset',get(gca,'TightInset'));
print(gcf,'localization_frame200_magnification','-dpng','-r300'); 
%%
[pixel_centers_row pixel_centers_col] = ndgrid(pixel_centers_y,pixel_centers_x);

localized_value_center = Localized(1:2);
localized_value_sigma_x = Localized(3);
localized_value_sigma_y = Localized(4);
localized_value_amplitude = Localized(5);
localized_value_background= Localized(6);


localized_image =1/(sqrt(localized_value_sigma_x*localized_value_sigma_y)*(2*pi))*localized_value_amplitude*exp(-(pixel_centers_col-...
    localized_value_center(1)).^2/(2*localized_value_sigma_x^2)).*exp(-(pixel_centers_row-...
    localized_value_center(2)).^2/(2*localized_value_sigma_y^2))+localized_value_background;


%%
figure()
subplot(121)
imagesc(ROI)
colormap gray
axis image
xlabel('Pixel ','Interpreter','Latex')
ylabel('Pixel ','Interpreter','Latex')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
subplot(122)
imagesc(localized_image)
colormap gray
axis image
xlabel('Pixel ','Interpreter','Latex')
ylabel('Pixel ','Interpreter','Latex')
set(gca,'YDir','reverse')
set(gca,'FontSize',18);
set(gca,'TickLabelInterpreter','latex');
ax=gca;
ax.LineWidth=1.5;
%print(gcf,'example_histogram','-dpng','-r300'); 

%%
save('frames_25000.mat');





