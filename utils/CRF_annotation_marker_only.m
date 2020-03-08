%%%% function to run CRF based annotation on marker only
%%%% Input examples
% data{1,33} = 'data_20190809_62g';
% data{1,34} = 'data_20190809_65g';
% data{1,35} = 'data_20190809_67g';
% numLandmarks = 0;
% numLabelRemove = 0;
% rng_fix_landmark = 42;
% num = 1;
% Note - don't forget to add path for specific CRF code

data = {};
data{1,1} = 'data_20190809_1g';
data{1,2} = 'data_20190809_3g';
data{1,3} = 'data_20190809_4g';
data{1,4} = 'data_20190809_5g';
data{1,5} = 'data_20190809_6g';
data{1,6} = 'data_20190809_7g';
data{1,7} = 'data_20190809_8g';
data{1,8} = 'data_20190809_9g';
data{1,9} = 'data_20190809_10g';
data{1,10} = 'data_20190809_11g';
data{1,11} = 'data_20190809_14g';
data{1,12} = 'data_20190809_19g';
data{1,13} = 'data_20190809_20g';
data{1,14} = 'data_20190809_22g';
data{1,15} = 'data_20190809_23g';
data{1,16} = 'data_20190809_24g';
data{1,17} = 'data_20190809_28g';
data{1,18} = 'data_20190809_29g';
data{1,19} = 'data_20190809_33g';
data{1,20} = 'data_20190809_34g';
data{1,21} = 'data_20190809_36g';
data{1,22} = 'data_20190809_37g';
data{1,23} = 'data_20190809_38g';
data{1,24} = 'data_20190809_39g';
data{1,25} = 'data_20190809_50g';
data{1,26} = 'data_20190809_52g';
data{1,27} = 'data_20190809_54g';
data{1,28} = 'data_20190809_55g';
data{1,29} = 'data_20190809_57g';
data{1,30} = 'data_20190809_58g';
data{1,31} = 'data_20190809_59g';
data{1,32} = 'data_20190809_60g';
data{1,33} = 'data_20190809_62g';
data{1,34} = 'data_20190809_65g';
data{1,35} = 'data_20190809_67g';
numLandmarks = 0;
numLabelRemove = 0;
rng_fix_landmark = 42;
num = 1;

addpath('C:\Users\Shivesh\Dropbox (GaTech)\PhD\GlobalBrainCode\GlobalBrainTrackingNew\functions\Annotation_CRF\NeuroPAL_validation\prediction_alternateStrains\AML5\RelPos')
match_score = struct();
for i = 1:size(data,2)
    [landmark_match_score,not_match,node_label,Neuron_head] = annotation_CRF_alternate_strain_3d_g_exsameLabel(data{1,i},numLandmarks,numLabelRemove,rng_fix_landmark,num);
    match_score(i).landmark_match_score = landmark_match_score;
    match_score(i).correct_fraction = max(landmark_match_score);
    match_score(i).not_match = not_match;
    match_score(i).node_label = node_label;
    match_score(i).Neuron_head = Neuron_head;
end