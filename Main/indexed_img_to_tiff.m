%%% function to convert indexed segmented image to tiff file so that it can
%%% be read in Vaa3D to check segmentation
%%% Can also be used to save any 3D image without providing cmap
function indexed_img_to_tiff(img,cmap,name)
    if ~isempty(cmap)
        imwrite(img(:,:,1),cmap,name)
        for i = 2:size(img,3)
            imwrite(img(:,:,i),cmap,name,'WriteMode','append')
        end
    else
        imwrite(img(:,:,1),name)
        for i = 2:size(img,3)
            imwrite(img(:,:,i),name,'WriteMode','append')
        end
    end
end