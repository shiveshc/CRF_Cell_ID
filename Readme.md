# CRF Cell ID

This software/tools were developed for automatically annotating biological identities of cells in dense brain images, visualizing automatic annotation results, and building new atlases using ground-truth annotations. In contrast to popular registration based methods, `CRF_Cell_ID` formulates a structured prediction problem that is solved with second-order Conditional Random Fields model. Thus `CRF_Cell_ID` searches for optimal labels for all cells that are most consistent with prior knowledge in terms of maintaining co-dependent features such as positional relationships. `CRF_Cell_ID` also provides a computationally scalable method for building atlases from annotated data.

To learn about details of the method, please read the <a href="https://www.biorxiv.org/content/10.1101/2020.03.10.986356v1">paper</a>. A detailed guide on using the framework is __coming soon__

This repository contains - 
1. Code for predicting identities in new datasets such as gene reporter expression, multi-cell calcium imaging or whole-brain imaging stacks
2. Code for building new atlas based on manually annotated datasets
3. Raw datasets used in the paper, ground-truth annotations in datasets and codes for reproducing results

<img src = "extra/readme_img_v2.jpg" width=50%>
	
	
	
  

  
