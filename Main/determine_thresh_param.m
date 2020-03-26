%%% function to determine thresh parameter for segmenting cells based on
%%% user input

function [thresh_param,accept_thresh,index_thresh_new] = determine_thresh_param(curr_img)
%     curr_img = imhmax(curr_img,0.01);
    loc_max = imregionalmax(curr_img);
    index = find(loc_max);
    [x,y,z] = ind2sub(size(loc_max),index);
    peaks = curr_img(index);
    thresh_param = input("Enter thresh_parameter value between 0-100. Higher values for detecting few but brightest cells - ");
    thresh_val = prctile(peaks,thresh_param);
    peak_thresh =  peaks(peaks>thresh_val);
    x_thresh = x(peaks>thresh_val);
    y_thresh = y(peaks>thresh_val);
    z_thresh = z(peaks>thresh_val); 
    index_thresh = index(peaks>thresh_val);
    [x_thresh_new,y_thresh_new,z_thresh_new,index_thresh_new] = remove_close_neurons(x_thresh,y_thresh,z_thresh,index_thresh,curr_img);
    figure,imshow(max(curr_img,[],3),[],'border','tight')
    caxis([0,0.6])
    hold on
    scatter(y_thresh,x_thresh,'.g')
    scatter(y_thresh_new,x_thresh_new,'.r')
    text(30,30,['Number of detected cells - ',num2str(size(x_thresh_new,1))],'Color','white')
    accept_thresh = input("Enter 'y' if accept thresh_parameter else enter 'n' - ");
    close all
end