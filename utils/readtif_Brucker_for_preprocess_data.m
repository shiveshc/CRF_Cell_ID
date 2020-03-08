%%% read brucker multichannel images
%%%%%%%%% Note %%%%%%%%%
%%% No need to supply one file (containing all timepoints) per channel
%%% If multiple cycles are present in this imaging session, be careful of
%%% what each cycle means. Current code works for if same number of z
%%% planes for each cycle

% addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\'))
num_channels = 1;
% in_direc = 'D:\PhD\Imaging\20190719_PS6058XOH10690\ZSeries-07192019-1400-219\';
% out_direc = 'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190719_PS6058XOH10690\1\';
% img_name = 'ZSeries-07192019-1400-219_Cycle00001_Ch2_000001.ome.tif'; %remember to change img_name
% vid_name = [];
% num_zplane = 20; % check from prarie reader or description file in the input folder
% divide = 1; % number of parts to divide image series to in case of huge time series (whole-brain video)

% in_direc = input_1{1,cnt};
% out_direc = input_2{1,n};
% img_name = input_3{1,cnt}; %remember to change img_name
% vid_name = [];
% num_zplane = input_4(1,n); % check from prarie reader or description file in the input folder
divide = 1; % number of parts to divide image series to in case of huge time series (whole-brain video)

%%% get number of cycles, number of images in each cycle and number of
%%% channels, img height and width
fileList = dir(in_direc);
cycle = {};
channel = {};
cnt = 1;
for i = 1:size(fileList,1)
    name = fileList(i).name;
    curr_cycle = regexp(name,'Cycle\d+','match');
    curr_channel = regexp(name,'Ch\d+','match');
    if ~isempty(curr_cycle) && ~isempty(curr_channel)
        cycle{cnt} = curr_cycle{1};
        channel{cnt} = curr_channel{1};
        cnt = cnt + 1;
    end
end
uniq_cycle = unique(cycle);
uniq_channel = unique(channel);
for i = 1:size(uniq_cycle,2)
    im_count(i) = sum(strcmp(cycle,uniq_cycle{i}));
    num_tp(i) = im_count(i)/(size(uniq_channel,2)*num_zplane);
end
[Height,Width] = size(imread([in_direc,img_name]));

%%% read stack images
if divide > 1
    for k = 1:divide
        for i = 1:size(uniq_cycle,2)
            for j = 1:size(uniq_channel,2)
                thisimage{j,i,k} = zeros(Height,Width,num_tp(i)/divide*num_zplane,'uint16');
            end
        end
    end
    
    counters = ones(size(thisimage));
    file_cnt = 0;
    curr_part = ones(size(uniq_channel,2),size(uniq_cycle,2));
    for i = 1:size(fileList,1)
        name = fileList(i).name
        curr_cycle = regexp(name,'Cycle\d+','match'); 
        curr_channel = regexp(name,'Ch\d+','match');
        if ~isempty(curr_cycle) && ~isempty(curr_channel)
            curr_cycle_index = find(strcmp(curr_cycle,uniq_cycle));
            curr_channel_index = find(strcmp(curr_channel,uniq_channel));
            curr_image = imread([in_direc,name]);
            if rem(counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)),num_tp(curr_cycle_index)/divide*num_zplane) == 0
                thisimage{curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)}(:,:,counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index))) = curr_image;
                counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)) = counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)) + 1;
                curr_part(curr_channel_index,curr_cycle_index) = curr_part(curr_channel_index,curr_cycle_index) + 1;
            else
                thisimage{curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)}(:,:,counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index))) = curr_image;
                counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)) = counters(curr_channel_index,curr_cycle_index,curr_part(curr_channel_index,curr_cycle_index)) + 1;
            end
        end
    end
    
    for k = 1:divide
        for i = 1:size(thisimage,1)
            for j = 1:size(thisimage,2) 
                thisimage{i,j,k} = reshape(thisimage{i,j},Height,Width,num_zplane,[]);
            end
        end
    end
else
    for i = 1:size(uniq_cycle,2)
        for j = 1:size(uniq_channel,2)
            thisimage{j,i} = zeros(Height,Width,im_count(i)/size(uniq_channel,2),'uint16');
        end
    end
    
    counters = ones(size(thisimage));
    for i = 1:size(fileList,1)
        name = fileList(i).name
        curr_cycle = regexp(name,'Cycle\d+','match'); 
        curr_channel = regexp(name,'Ch\d+','match');
        if ~isempty(curr_cycle) && ~isempty(curr_channel)
            curr_cycle_index = find(strcmp(curr_cycle,uniq_cycle));
            curr_channel_index = find(strcmp(curr_channel,uniq_channel));
            curr_image = imread([in_direc,name]);
            thisimage{curr_channel_index,curr_cycle_index}(:,:,counters(curr_channel_index,curr_cycle_index)) = curr_image;
%             thisimage{curr_channel_index,curr_cycle_index}(:,:,num_zplane - counters(curr_channel_index,curr_cycle_index) + 1) = curr_image;
            counters(curr_channel_index,curr_cycle_index) = counters(curr_channel_index,curr_cycle_index) + 1;
        end
    end
    
    for i = 1:size(thisimage,1)
        for j = 1:size(thisimage,2) 
            thisimage{i,j} = reshape(thisimage{i,j},Height,Width,num_zplane,[]);
        end
    end
end

%%% compile final images consisting of all cycles. Be careful of what each
%%% cycle means
finalimage = cell(1,size(uniq_channel,2));
for i = 1:size(uniq_channel,2)
    for j = 1:1
        finalimage{i} = thisimage{i,j};
    end
end
for i = 1:size(uniq_channel,2)
    for j = 2:size(uniq_cycle,2)
        finalimage{i} = cat(4,finalimage{i},thisimage{i,j});
    end
end

%%% optional - save individual channel stacks compatible with Vaa3D
% indexed_img_to_tiff(thisimage{1,1},[],strcat(out_direc,'thisimage_g.tif'))
%     % indexed_img_to_tiff(thisimage_r,[],strcat(out_direc,'\Vaa3d_',file_r))
% 
% %%% optional - save channels for neuron activity GUI
% % name_g = [out_direc,'RFP.tif'];
% % thisimage_g = finalimage{1,1};
% % imwrite(max(thisimage_g(:,:,:,1),[],3),name_g)
% % for k = 2:size(thisimage_g,4)
% %     imwrite(max(thisimage_g(:,:,:,k),[],3),name_g,'WriteMode','append');
% % end
% 
% %%% create video (consists of all cycles)
% g_caxis_scaling = [0,0.3];
% r_caxis_scaling = [0,0.3];
% 
% writerObj = VideoWriter([out_direc,vid_name]); % Name it.
% writerObj.FrameRate = 25; % How many frames per second.
% open(writerObj);
% cnt = 1;
% if num_channels > 1
%     if divide > 1
%         for i = 1:size(uniq_cycle,2)
%             for k = 1:divide
%                 for j = 1:size(thisimage{1,i,k},4)
%                     fig_g = mat2gray(max(thisimage{1,i,k}(:,:,:,j),[],3));
%                     fig_r = mat2gray(max(thisimage{2,i,k}(:,:,:,j),[],3));
%         %             fig_g_norm = (fig_g - min_p(1))/(max_p(1) - min_p(1));
%                     fig_g_norm_adj = imadjust(fig_g,[g_caxis_scaling(1,1);g_caxis_scaling(1,2)],[0;1]);
%         %             fig_r_norm = (fig_r - min_p(1))/(max_p(1) - min_p(1));
%                     fig_r_norm_adj = imadjust(fig_r,[r_caxis_scaling(1,1);r_caxis_scaling(1,2)],[0;1]);
%                     fig_b_norm = zeros(size(fig_g));
%                     rgb_img = cat(3,fig_r_norm_adj,fig_g_norm_adj,fig_b_norm);
%                     rgb_img(30:31,330:350,:) = 1;               % scale bar for 5 um (0.23 um/pixel for opterra)
% 
%                     imshow(rgb_img)
%                     text(10,10,['t = ',num2str(cnt),' cycle = ',num2str(i), ' part = ',num2str(k)],'Color','white')
%         %             print(gcf,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode_Charles\GlobalBrainTrackingNew\20171006_MEC_Device_AML32_Brucker\AML32_3\Frames\',num2str(i)],'-dpng','-r300')
%                     frame = getframe(gcf);
%                     writeVideo(writerObj, frame);
%                     cnt = cnt + 1;
%                 end
%             end
%         end
%     else
%         for i = 1:size(uniq_cycle,2)
%             for j = 1:size(thisimage{1,i},4)
%                 fig_g = mat2gray(max(thisimage{1,i}(:,:,:,j),[],3));
%                 fig_r = mat2gray(max(thisimage{2,i}(:,:,:,j),[],3));
%     %             fig_g_norm = (fig_g - min_p(1))/(max_p(1) - min_p(1));
%                 fig_g_norm_adj = imadjust(fig_g,[g_caxis_scaling(1,1);g_caxis_scaling(1,2)],[0;1]);
%     %             fig_r_norm = (fig_r - min_p(1))/(max_p(1) - min_p(1));
%                 fig_r_norm_adj = imadjust(fig_r,[r_caxis_scaling(1,1);r_caxis_scaling(1,2)],[0;1]);
%                 fig_b_norm = zeros(size(fig_g));
%                 rgb_img = cat(3,fig_r_norm_adj,fig_g_norm_adj,fig_b_norm);
%                 rgb_img(30:31,330:350,:) = 1;               % scale bar for 5 um (0.23 um/pixel for opterra)
% 
%                 imshow(rgb_img)
%                 text(10,10,['t = ',num2str(cnt),' cycle = ',num2str(i)],'Color','white')
%     %             print(gcf,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode_Charles\GlobalBrainTrackingNew\20171006_MEC_Device_AML32_Brucker\AML32_3\Frames\',num2str(i)],'-dpng','-r300')
%                 frame = getframe(gcf);
%                 writeVideo(writerObj, frame);
%                 cnt = cnt + 1;
%             end
%         end
% 
%     end
% else
%     if divide > 1
%         for i = 1:size(uniq_cycle,2)
%             for k = 1:divide
%                 for j = 1:size(thisimage{1,i,k},4)
%                     fig_g = mat2gray(max(thisimage{1,i,k}(:,:,:,j),[],3));
%         %             fig_g_norm = (fig_g - min_p(1))/(max_p(1) - min_p(1));
%                     fig_g_norm_adj = imadjust(fig_g,[g_caxis_scaling(1,1);g_caxis_scaling(1,2)],[0;1]);
%                     fig_r_norm_adj = zeros(size(fig_g));
%                     fig_b_norm = zeros(size(fig_g));
%                     rgb_img = cat(3,fig_r_norm_adj,fig_g_norm_adj,fig_b_norm);
%                     rgb_img(30:31,330:350,:) = 1;               % scale bar for 5 um (0.23 um/pixel for opterra)
% 
%                     imshow(rgb_img)
%                     text(10,10,['t = ',num2str(cnt),' cycle = ',num2str(i), ' part = ',num2str(k)],'Color','white')
%         %             print(gcf,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode_Charles\GlobalBrainTrackingNew\20171006_MEC_Device_AML32_Brucker\AML32_3\Frames\',num2str(i)],'-dpng','-r300')
%                     frame = getframe(gcf);
%                     writeVideo(writerObj, frame);
%                     cnt = cnt + 1;
%                 end
%             end
%         end
%     else
%         for i = 1:size(uniq_cycle,2)
%             for j = 1:size(thisimage{1,i},4)
%                 fig_g = mat2gray(max(thisimage{1,i}(:,:,:,j),[],3));
%     %             fig_g_norm = (fig_g - min_p(1))/(max_p(1) - min_p(1));
%                 fig_g_norm_adj = imadjust(fig_g,[g_caxis_scaling(1,1);g_caxis_scaling(1,2)],[0;1]);
%                 fig_r_norm_adj = zeros(size(fig_g));
%                 fig_b_norm = zeros(size(fig_g));
%                 rgb_img = cat(3,fig_r_norm_adj,fig_g_norm_adj,fig_b_norm);
%                 rgb_img(30:31,330:350,:) = 1;               % scale bar for 5 um (0.23 um/pixel for opterra)
% 
%                 imshow(rgb_img)
%                 text(10,10,['t = ',num2str(cnt),' cycle = ',num2str(i)],'Color','white')
%     %             print(gcf,['C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode_Charles\GlobalBrainTrackingNew\20171006_MEC_Device_AML32_Brucker\AML32_3\Frames\',num2str(i)],'-dpng','-r300')
%                 frame = getframe(gcf);
%                 writeVideo(writerObj, frame);
%                 cnt = cnt + 1;
%             end
%         end
%     end
% end
% close(writerObj)