# The Influence of Nuisance Parameter Uncertainty on Parameter Estimation in Practical Data Science Models

* Abstract

The uncertainty associated with nuisance parameters profoundly influences the estimation of main parameters, thereby influencing predictive outcomes. This paper delves into the exploration of this uncertainty in the context of univariate nuisance parameter GARCH models and multivariate nuisance parameter Neural Network models. Leveraging Spall's theoretical framework for nuisance parameter uncertainty, we construct practical models that comprehensively capture the variance of primary parameters, considering the inherent uncertainty introduced by nuisance parameters. Our investigation reveals that the impact of nuisance parameters' uncertainty varies across different model structures, even within the same model. Furthermore, our study affirms the practical effectiveness of the adjusted models, which not only provide more credible confidence ellipses and intervals but also enhance the reliability of predictive analytics.

* Introduction

Mathematical formulas serve as the foundation for contemporary models across various disciplines. The process of model establishment fundamentally involves parameter estimation. To optimize computational efficiency in terms of both calculation quantity and time, researchers categorize parameters that demand substantial time for calculation within the model but are not the focus of the study as nuisance parameters. These parameters are then assigned values based on empirical knowledge to streamline the calculation process. Despite being considered nuisances, these parameters are integral to the model, and their uncertainties directly influence the estimation of primary parameters, consequently impacting the model's accuracy.

The properties of primary parameters in the context of uncertainty surrounding assumed values for nuisance parameters was scrutinized by Spall(1989). The formula introduced provides a more accurate representation of the uncertainty associated with nuisance parameters in the variance of primary parameters. This enhancement enables confidence intervals to encompass the true value of the primary parameter.

Building upon Spall's theory, this article applies it to two distinct data science models: 1) the GARCH model and 2) Neural Network, situated in the realms of finance and machine learning, respectively. The objective is to investigate how the uncertainty surrounding nuisance parameters affects these models individually.
  
