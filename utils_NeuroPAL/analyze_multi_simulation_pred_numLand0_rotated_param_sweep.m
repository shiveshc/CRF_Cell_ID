%%%% function to analyze probabilistic label run of neuropal data by
%%%% removing random combination of labels in each run. This function is
%%%% for 0 number of landmarks. Use 'analyze_multi_simulation_pred.m' in
%%%% case of landmarks. This function only for rotated strains

%%%% Input - input directories of probabilistic label results

function pred_struct = analyze_multi_simulation_pred_numLand0_rotated_param_sweep()
%     addpath(genpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\'))
    in_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190705_27\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190706_16\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190706_20\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190710_21\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_1\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_4\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_8\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_9\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_10\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_11\',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\loo_posboth_colrot\ensem_int_colconst\nonpca\multi_data_20190720_13\'};
    
    data = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190705_27_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190706_16_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190706_20_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190710_21_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_1_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_4_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_8_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_9_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_10_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_11_pca_ensem_int_colconstnorm.mat',...
    'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_NeuroPAL\RotatedStrains\RelPos_Col\data_20190720_13_pca_ensem_int_colconstnorm.mat'};

    marker_direc = {'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190705_OH15495_array\27',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\16',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190706_OH15495_array\20',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190710_OH15495_array\21',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\1',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\4',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\8',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\9',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\10',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\11',...
        'C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\20190720_OH15495_array\13'};
    
    pred_struct = struct();
    param_sweep = cat(2,[0.1:0.1:0.9]',[0.9:-0.1:0.1]');
    
    for f = 1:size(in_direc,2)
        top_correct_frac = [];
        top_3_correct_frac = [];
        top_5_correct_frac = [];
        fileList = dir(in_direc{1,f});
%         load(data{1,f})
        [X,Y,Z,marker_name,marker_index] = read_marker_files(marker_direc{1,f});
        landmark_list = marker_index;
%         landmark_list = [1:1:size(marker_index,1)]';
        
        for p = 1:size(param_sweep,1)
            labels = {};
            cnt = 1;
            for i = 1:size(fileList,1)
                name = fileList(i).name;
                if ~strcmp(name,'.') && ~strcmp(name,'..')
                    underscore_pos = strfind(name,'_');
                    param1 = str2num(name(1,underscore_pos(1,end-2)+1:underscore_pos(1,end-1)-1));
                    param2 = str2num(name(1,underscore_pos(1,end-1)+1:underscore_pos(1,end)-1));
                    if param1 == round(param_sweep(p,1),1)
                        load([in_direc{1,f},name])

                        landmark_to_be_pred = landmark_list;
                        for k = 1:size(landmark_to_be_pred,1)
                            if experiments(1).node_label(landmark_to_be_pred(k,1),1) == 0
                                labels{k,cnt} = '';
                            else
                                labels{k,cnt} = experiments(1).Neuron_head{experiments(1).node_label(landmark_to_be_pred(k,1),1),1};
                            end
                        end
                        cnt = cnt + 1;
                    end
                end
            end
            if ~isempty(labels)
                top_labels = calculate_top_labels(labels);
                true_labels = marker_name;
                top_correct = 0;
                top_3_correct = 0;
                top_5_correct = 0;
                for tr_lb = 1:size(true_labels,1)
                    if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(tr_lb).top_labels(1,1))))
                        top_correct = top_correct + 1;
                    end
                    if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(tr_lb).top_labels(1,1:min(3,size(top_labels(tr_lb).top_labels,2))))))
                        top_3_correct = top_3_correct + 1;
                    end
                    if ~isempty(find(strcmp(true_labels(tr_lb,1),top_labels(tr_lb).top_labels(1,1:min(5,size(top_labels(tr_lb).top_labels,2))))))
                        top_5_correct = top_5_correct + 1;
                    end
                end
                top_correct_frac = [top_correct_frac;param_sweep(p,1),param_sweep(p,2),top_correct/size(top_labels,2)];
                top_3_correct_frac = [top_3_correct_frac;param_sweep(p,1),param_sweep(p,2),top_3_correct/size(top_labels,2)];
                top_5_correct_frac = [top_5_correct_frac;param_sweep(p,1),param_sweep(p,2),top_5_correct/size(top_labels,2)];
            else
                top_correct_frac = [top_correct_frac;param_sweep(p,1),param_sweep(p,2),-1];
                top_3_correct_frac = [top_3_correct_frac;param_sweep(p,1),param_sweep(p,2),-1];
                top_5_correct_frac = [top_5_correct_frac;param_sweep(p,1),param_sweep(p,2),-1];
            end
        end

        pred_struct(f).data_name = data{1,f};
        pred_struct(f).top_correct_frac = top_correct_frac;
        pred_struct(f).top_3_correct_frac = top_3_correct_frac;
        pred_struct(f).top_5_correct_frac = top_5_correct_frac;
    end
end

function top_labels = calculate_top_labels(labels)
    top_labels = struct();
    for i = 1:size(labels,1)
        uniq_labels = unique(labels(i,:));
        count = [];
        for j = 1:size(uniq_labels,2)
            count(1,j) = sum(strcmp(uniq_labels(1,j),labels(i,:)));
        end
        [sort_count,sort_index] = sort(count,'descend');
        top_labels(i).top_labels = uniq_labels(1,sort_index(1,1:min(5,size(sort_index,2))));
    end
end