# Monte Carlo Simulation Analysis of Fama–French Five-Factor Model

## Abstraction

Linear Regression with O.L.S method is most popular approach for modelling the relationship between variables in academic and practical fields. However, researchers and practitioners often misuse the O.L.S properties by ignoring the prerequisites. In this paper, we use Monte Carlo simulations to validate the Fama–French five-factor model and the O.L.S properties. We confirm the validity of five factors in the Fama–French model, and uncover six interesting theoretical findings of O.L.S properties.

## Conclusion

In this study, we firstly do O.L.S to investigate the factors to describe stock returns (i.e. five factors in the Fama– French model). We downloaded data from Jan 1st, 2008 to Dec 29th, 2017 from the Kenneth R. French homepage. Then we got the true corresponding coefficients of each factors by O.L.S regression using real data. After that, we utilize Monte Carlo simulations to test the validity of the Fama–French five-factor model and prove some properties of OLS estimator (such as unbiasedness), confidence interval, t-test, Type I error and model selection.

We confirm the validity of five factors in the Fama–French model. In addition, we compare the relative importance between these factors, and find that CAPM is indeed the most important explanatory, and HML after it. These findings provide a good deal of insight into portfolio management and asset pricing.

Our results also uncovered six interesting theoretical findings through the regression model, 

1. The O.L.S. estimators of the unknown coefficients and error variance are unbiased. 
2. The “correct” meaning of a $100(1-\alpha)$% confidence interval of an unknown coefficient is that roughly the true coefficient value of the coefficient will fall into the generated confidence interval in $100(1-\alpha)$% replications.
3. The significance level of the t test for testing a linear hypothesis concerning one or more coefficients is the probability of committing a Type I error.
4. The t test is unbiased.
5. if some of the important variables are omitted, the estimators of the remaining coefficients will become biased.
6. in model selection, practitioners should keep balance between unbiased estimator and low variance of estimator by using some model selection techniques.

In conclusion, the Monte Carlo simulations help us validate each factor in the Fama–French five-factor model and understand the properties of O.L.S. In the future, we will use O.L.S regression to explore more factors in other asset pricing model.

**Note**: This is the code for project of FB8916@CityU, Semester B, 2017-2018. The topic is about using Monte Carlo technique to investigate some linear regression properties of the Fama-French five-factor model.


Details can be found on post [Monte Carlo Experiment for Fama French 5 Factor Model](http://blog.baoduge.com/MC_FF5/).
