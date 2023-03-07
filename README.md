# MLE_Continuation
Matlab code to reproduce the figures in the preprint: 

The code is organized into three folders corresponding to the Phenotype PDE and Viral dynamics examples as well as generic template. In each example folder, there is a main driver code use to reproduce the corresponding figure from the manuscript. In addition, there are two helper functions. These helper functions are the log-likelihoood, or objective function, that is minimized during MLE and a script to output model values at the sample points t_i for given parameters \theta. 


To reproduce Figure 1, run:

\MLE_Continuation\Phenotype PDE\Continuation Technique\PLOS_CompBio_Predictor_Continuation_4Param_Treated_Figure 

To reproduce Figure 2, run:

\MLE_Continuation\Phenotype PDE\Experimental Prediction\PLOS_CompBio_ExperimentalPrediction_Example  

To reproduce Figure 3, run:

\MLE_Continuation\Viral dynamics\HIV_Example

