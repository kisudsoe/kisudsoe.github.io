<a href="https://kisudsoe.github.io"><img src="img/favicon.png" width="30px" /></a>[Home](https://kisudsoe.github.io)

`2017-03-25` Posted | `2018-06-24` Add description

[TOC]

# Confounding factor control

Here, we are going to watch how to control confounding factors.



## 1. Stratification

In your data, **bias** make things complicate. To avoid the bias, we can try three different ways such as sampling randomization, restricted sampling and matched sampling.

### Simpson's paradox

The **danger of bias** can be displayed from **"Simpson's paradox"**, which is also known as **reversal paradox**. This example shows that we would misunderstand as opposite trend against the truth.

![](/img/2017-03-25-Confounding_factor_control/슬라이드2.JPG)



### Stratification

Your data could be grouped by some ways, which we called as **stratification**. Using this way appropriately, we could find difference between sub-groups in your data. Again, you can find bias using this way!


![](/img/2017-03-25-Confounding_factor_control/슬라이드3.JPG)



However, **stratification** have limitation from sample number. When you are dividing data into many sub-groups, each sub-group has only few sample number, which make loss of statistical power.

![](/img/2017-03-25-Confounding_factor_control/슬라이드4.JPG)



## 2. Multivariate models

Another way to avoid confounding factor is using multivariate models, including logistic regression, linear regression and analysis of covariance (ANCOVA). Among them, I prefer using linear regression model which is intuitive to understand and good for visualization as plots.

![](/img/2017-03-25-Confounding_factor_control/슬라이드5.JPG)



### Finding fitted model

To generate linear regression model, you can add or eliminate variables. By some changes, you can find the most fitted model. To find the most fitted model, you can consider the probability (p) of model. In fact, lower p-value means more fitted.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드6.JPG)



### Regression model types

Regression model has various types such as univariate and multivariate by number of outcome variable, simple and multiple by number of predictor variable as well as linear and nonlinear by mathematical formula. In addition, analysis of variance (ANOVA), ANCOVA and Logistic regression are also include regression model.

When you have trouble to find the best fitted model, you should try to generate different types of models as well as control different variables by comparing the p-values of each model.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드7.JPG)



### Example

![2](/img/2017-03-25-Confounding_factor_control/슬라이드8.JPG)



## 3. Linear model selection

You might have curiosity to interaction of variables. Again, to identify their fitness for your data, you can compare the p-values between models with or without interaction of variables. Another way to compare model fitness is using **likelihood ratio test** between models.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드9.JPG)



### Likelihood ratio test

Here, I show you a R package to perform likelihood ratio test, which is lrtest function in 'epiDisplay'.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드10.JPG)



### Choose explanatory variables

There are three different methods to generate models by controlling variables. **Backward elimination** method start full model, which use every variables for model. Then eliminate one by one variable from the model. **Forward selection** method is opposite way of the Backward elimination. And **stepwise method** is started from constant number only model. And then you add a variable to model and compare other variable models. And then you choose a model with the smallest p-value. Then add another variable to the model and do the same things again.

Whatever you can choose a method among the three. But you have to invest enough time to find the best fitted model to your data.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드11.JPG)



### Example: Results

![2](/img/2017-03-25-Confounding_factor_control/슬라이드13.JPG)



## 4. Heritability (*h*^2^)

Calculation of heritability is another story of regression model. I just add this section to refer it later.

![2](/img/2017-03-25-Confounding_factor_control/슬라이드14.JPG)



### Example: Results

Here is an example for the linear regression by Hoffman JM *et al.* (2014).

![2](/img/2017-03-25-Confounding_factor_control/슬라이드15.JPG)



![2](/img/2017-03-25-Confounding_factor_control/슬라이드16.JPG)



![2](/img/2017-03-25-Confounding_factor_control/슬라이드17.JPG)

---

<div id="disqus_thread"></div>
<script>

(function() { // DON'T EDIT BELOW THIS LINE
var d = document, s = d.createElement('script');
s.src = 'https://kisudsoe-github-io.disqus.com/embed.js';
s.setAttribute('data-timestamp', +new Date());
(d.head || d.body).appendChild(s);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript">comments powered by Disqus.</a></noscript>
<script id="dsq-count-scr" src="//kisudsoe-github-io.disqus.com/count.js" async></script>