# CICCplot

CICCplot is an R function that can be used to plot Conditional Item Characteristic Curves. These are useful when fitting a Rasch model to datainnorder to graphically detect misfit of an item in the model. A CICC describies the expected outcome of an item when conditioning on the total score of all items in the model. When plotting CICC's the total score will be on the x-axis and the conditional expected item score will be on the y-axis. For dichotomous items the expected item score corresponds to the probability of getting item score '1' and is always zero when the total score is zero and one when the total score is the maximal possible score. Between these to points the curve will be monotone increasing. 

## Before using CICCplot

Before the `CICCplot`-function can be used, the following packages need to be installed (if not already done) and included by using `library()`.

```
library(eRm)
library(iarm)
library(ggplot2)
library(ggpubr)
```

to use the function `CICCplot` directly from a R-script containig the function. Therefore it is necessary to source the function into our current work, when it should be applied.

```
source("CICCplot_source.R")
```

## Plotting CICC's in R

How the function `CICCplot` can be used to compare model CICC's to observed data will now be illustrated by two data examples.

### AMTS

Abbreviated Mental Test Score (AMTS) is a score used to detect dementia. It consists of 10 questions (taking values 0 and 1), which can either be correctly answered (1 point) or incorrectly answered (0 points).  In the `iarm` package a dataset called `amts` containing responses to AMTS from 197 patients is available. Beside the scores of the 10 questions, the data also contain an identification number, agegroup and sex of each patient. In the following analysis we are only interested in the answers to the 10 AMTS-questions, and the remainig variables will be left out of the analysis.

```
library(iarm)
head(amts)
```

First we fit a dichotomous Rasch model to the data using `RM()`-function from the `eRm` package. When fitting the model, we include the answers to all the 10 questions in AMTS. 

```
mod1 <- RM(amts[,4:13], sum0 = FALSE)
```

The estimated parameters can be found using `summary()`.

```
summary(mod1)
```

`CICCplot` can now be used to investigate misfit of the model to the data for each item. Different arguments can be stated in the function to create the plot that the user needs. The following three arguments always needs to be specified:

- `model` : A model object passed from the RM()-function in the eRm-package.
- `which.item` : An integer or vector giving the item(s), for which a CICC-plot should be constructed. Default is `which.item = 1`. The argumet will not be used if `all.items = TRUE`.
- `lower.groups` : A vector used for dividing the set of possible total scores into intervals, for which the emperical expected item-score will be calculated and added to the plot. The vector should contain the lower points of the intervals, that the set of possible total scores should be divided into. If zero doesn't appear in the vector, it will be added automaticly. If `lower.groups = "all" `, the emperical expected item-score will plottet for every possible total score. 

Constructing a plot for item 2 in `mod1`, where the observed mean of item scores is included for every possible total score can be done as follows:

```
CICCplot(mod1, which.item = 2, lower.groups = "all")
```

It can be useful to be able to divide the total scores into groups, for which the observed outcomes will then be plotted. This can for instance be if the number of observations is small for some total scores. For the `amts` data we calculate the total score for each individual by `rowSums()` and the find the total counts for each total score by using `table()`.

