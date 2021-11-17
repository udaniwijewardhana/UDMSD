# UDMSD

## Single species annual spatial and temporal joint model to geostatistical zero-inflated data

To estimate the persistence, user can consider INLA models with zero-inflation probability model which means joint models. Joint temporal, spatial and spatio-temporal models are especially developed for data which has access zeros. Access zeros are divided into two types such as structural zeros and random zeros. Structural zeros refer to zero responses by those subjects whose count response will always be zero and random zeros that occur to subjects whose count response can be greater than zero but appear to be zero due to sampling variability (He et al., 2014). The responses of species data came from two different distributions such as occurrence and abundance which have two models for each response that are affected by spatial and temporal common factors. Therefore, it is better to use joint model. Negative Binomial model is more in line for zero inflation species data which resolve the overdispersion issue as well. Therefore in this window users can get the summary outputs for Joint NB model (Cameron and Trivedi, 2013), Joint Hurdle NB model (Cragg, 1971; Mullahy, 1986) or Joint Zero Inflated Negative Binomial model (Cameron and Trivedi, 1998) which are the most common zero inflation joint models to identify the significance of the predictors and to identify the persistency.

A Bayesian hierarchical modelling approach is used here to allow us to conveniently account for structures in parameter uncertainty and potential dependence such as spatial and temporal structures. With a Bayesian approach, a joint posterior distribution is obtained for the process and parameters of the given data. This is conducted using Integrated Nested Laplace Approximation (INLA). INLA is popular as an approximation tool for fitting Bayesian models. INLA is an alternative robust method for the traditional Markov Chain Monte Carlo (MCMC) based Bayesian analyses (Paul et al. 2010). INLA approach is a numerically implemented analytical solution for approximating posterior marginals in hierarchical models with latent Gaussian processes. The key advantages of INLA are the ease with which complex models can be created and modified, without the need to write complex code, and the speed at which inference can be done even for spatial problems with hundreds of thousands of observations (Sarul, 2015). In spatio-temporal settings, it is often assumed that the covariance is separable in space and time, and thus, the temporal structure may be modelled using an auto-regressive process (Cressie, 2011). Fitted INLA with log-linear regression for the positive counts (truncated count model) and a logistic regression for the zero counts using a Bernoulli distribution in spatial and temporal scales. When we are dealing with process defined over a continuous domain, one can express a large class of random fields as solution of stochastic partial differential equations (SPDEs). In R-INLA this solution is approximated using high dimensional basis representation with simple local basis function. 

These basis functions are defined over a triangulation of the domain; this triangulation is the mesh https://ourcodingclub.github.io/tutorials/inla/). This app can only create non-convex meshes. The general hierarchical model can be described in three phases such as data model, process model and parameter model. Adapted these INLA models using R-INLA (Martins, Simpson, Lindgren and Rue, 2013). 

To know if the rate at which abundance is changing over time differs according to a relevant predictor variable, we have included the facility to add interaction terms between any two predictor variables in our regression models. The interaction with a categorical variable tells us what the difference in slope is and whether this difference is significant.

UDMSD is a Shiny web application that allows to visualize spatial, temporal or spatio-temporal abundance/occurrence data, estimate significant numerical factors and the trend. It is addressed to ecologists interested in analysing species abundance/occurrence data but lacking the appropriate theoretical knowledge to use this statistical software is required. The application allows to fit Bayesian Zero inflation models to obtain significant predictor estimates and their uncertainty by using INLA. To carry out these analyses users simply need to click the buttons that create the input files required, execute the software and process the output to generate tables of values and plots with the results. The application allows user interaction and creates interactive visualizations such as species distribution map and posterior maps. This INLA models could be done only for geostatistical data. 

The application consists of 2 pages with main window: main window allows the user to upload the input files (data file and predictors file) and also gives the option to normalize or standardize the data or predictors data. 

1) Data display page where the user can see the annual data table of the uploaded data file and uploaded annual predictors table.

2) an analysis page where statistical analyses are carried out, and summary results and posterior plots can be visualized.

In this page, user can fit joint species distribution models using INLA. Since species abundance data usually have access and true zeros, zero inflation models are ideal to fit the abundance. Therefore, while overcoming overdispersion Negative Binomial model, Negative Binomial Hurdle model and ZINB model have given as options using log linear regression as the link function. Binomial distribution has used for occurrence using logistic regression as the link function. First, a user must adjust the parameters to build the INLA mesh. Only non-convex meshes could be built in this app. Then users can fit spatial, temporal or spatio-temporal joint species distribution while changing the model options and distribution. Users also have the option to select which type of spatial and temporal models could be done. Temporal model has the option to choose the random effect model such as 'iid', 'ar1', 'rw1' or 'rw2' and the spatial effect model such as 'iid' for temporal models and 'spde' for spatial or spatio-temporal models. Posterior plots are visualizing for spatial or spatio-temporal models.

# Input file: 

Dataframe with "count", "zero", Year", "Latitude", "Longitude", "Locality" and other numeric predictors (optional). These variables are case sensitive. Only applicable up to five predictor variables. Predictor variables display as "p.z1", "p.y2", etc. The order is the exact order in dataframe.    

1. "Locality" - Detected Location 
2. "Latitude" - Latitude coordinate of the Location 
3. "Longitude" - Longitude coordinate of the Location 
4. "Year" - Detected Year
5. count - Dependent variable for count model
6. zero - Dependent variable for zero inflation model

The 'Single-species Joint Model' gives a combination of count model and zero-inflation model results with spatial and temporal scale using INLA. This window also visualizes posterior mean and standard deviation plots for spatial and spatio-temporal models. This shows INLA mesh, summary results of the model and posterior plots (when applicable). Here, users can change the mesh parameters as well as spatial and temporal effect model in order to find the most suitable model. 

### Installation

To build this Shiny app, we need to clone the Zip file from UDMSD and save it in our computer. This folder contains a sample data .CSV file, the vignette and app.R file. Then, we can launch the app by clicking the Run App button at the top of the RStudio editor or by executing runApp("appdir_path")where appdir_path is the path of the directory that contains the app.R file. For this we need to install R and RStudio in our computer. The users who do not have R in their computer can use UDMSD to launch the Shiny app. The application uses the R-INLA package which can be downloaded from http://www.r-inla.org/download. 

```r
shiny::runGitHub( "UDMSD", "uwijewardhana") 
```