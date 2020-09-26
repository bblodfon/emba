---
title: 'emba: R package for analysis and visualization of biomarkers in boolean model ensembles'
tags:
  - R
  - boolean networks
  - logical modeling
  - biomarkers
  - mechanistic models
  - drug synergies
  - anti-cancer drug combinations
  - druglogics
authors:
  - name: John Zobolas
    orcid: 0000-0002-3609-8674
    affiliation: 1, 2
  - name: Martin Kuiper
    orcid: 0000-0002-1171-9876
    affiliation: 1
  - name: Åsmund Flobak
    orcid: 0000-0002-3357-425X
    affiliation: 2, 3
affiliations:
 - name: Department of Biology, Norwegian University of Science and Technology (NTNU), Trondheim, Norway
   index: 1
 - name: Department of Clinical and Molecular Medicine, Norwegian University of Science and Technology (NTNU), Trondheim, Norway
   index: 2
 - name: The Cancer Clinic, St. Olav’s Hospital, Trondheim, Norway
   index: 3
date: 31 July 2020
bibliography: paper.bib
---

# Introduction

Computational modeling of cellular systems has been one of the most powerful tools used to build interpretable knowledge of biological processes and help identify molecular mechanisms that drive diseases such as cancer [@Aldridge2006].
In particular, the use of logical modeling has proven to be a substantially useful approach, since it allows the easy construction, simulation and analysis of predictive models, capable of providing a qualitative and insightful view on the extremely complex landscape of biological systems [@Morris2010; @Wang2012; @Abou-Jaoude2016].
Such mechanistic models, with the systematic integration of prior knowledge and experimental data, have been extensively used to better understand what drives deregulation of signal transduction, the outcome of which is the manifestation of diseases [@Traynard2017].
Furthermore, their explanatory power has been used to provide insights into a drug’s mode of action, investigate the mechanisms of resistance to drugs [@Eduati2017] and suggest new therapeutic combination candidates, among others [@Flobak2015].

One of the major challenges in systems medicine, has been the identification of scientifically validated, predictive biomarkers that correlate with patient response to given therapies.
The analysis of biological predictive markers of pharmacologic response can not only further our understanding of the systemic processes involved in diseases but can also help to classify patients into groups with similar responses to specific therapeutic interventions, advancing personalized medicine [@Senft2017].
In addition, the identification of biomarkers in tumor cells (e.g. mutations) has enabled the discovery of drug targets which are utilized in combinatorial molecular-targeted therapies - a strategy which aims to treat specific patient subgroups and has shown larger overall survival rates and reduced side-effects than monotherapy [@Al-Lazikani2012].
Despite the huge advancements towards drug combination therapy, genetic heterogeneity, drug resistance and drug combination synergy mechanisms still pose fundamental challenges to clinicians, modelers and lab researchers.

To help bridge the model simulation results with the (clinical) laboratory observations, several optimization methods have been used, such as model calibration, parameter estimation and sensitivity analysis.
These methods also allow us to determine which model parameters have the biggest influence in the overall behaviour of the system [@Aldridge2006].
For example, in @Frohlich2018, a computational framework that allowed for the efficient parameterization and contextualization of a large-scale cancer signaling network, was used to predict combination treatment outcome from single drug data.
This model was calibrated to fit and accurately describe specific cell-line experimental data, while enabling the identification of biomarkers of drug sensitivity as well as molecular mechanisms that affect drug resistance.
Furthermore, in @Dorier2016, a network optimization approach which topologically parameterized boolean models according to a genetic algorithm was used, in order to best match the experimentally observed behaviour.
This method resulted in an ensemble of boolean models which can be used to simulate response under drug perturbations in order to assess the underlying mechanisms and to generate new testable hypotheses.
Such an aggregation of best-fit models (wisdom of the crowds) has been shown to be quite robust and effective for model prediction performance [@Marbach2012].

# Statement of need

There is a plethora of software tools devoted to the qualitative modeling and analysis of biological networks.
The Consortium for the development of Logical Models and Tools (CoLoMoTo) is a community effort which aims to standardize the representation of logical networks and provide a common repository of methods and tools to analyze these networks [@Naldi2015].
Furthermore, to facilitate the access to several software logical modeling tools and enable reproducible computational workflows, the CoLoMoTo Interactive Notebook was introduced as a unified computational framework [@Naldi2018a].
The incorporated tools are accessed via a common programming interface (though originally implemented in different programming languages e.g. Java, Python, C\texttt{++} and R) and offer a collection of features like accessing online model repositories [@Helikar2012], model editing [@Naldi2018b], dynamical analysis (finding attractors, stochastic simulations, reachability properties, model-checking techniques) [@Mussel2010; @Klarner2016; @Stoll2017; @Pauleve2017; @Naldi2018] and model parameterization/optimization to fit perturbation signaling data [@Terfve2012; @Gjerga2020].
Despite the diverse and multi-purpose logical modeling tools that exist, there is still a lack of data analysis-oriented software that assists with the discovery of predictive biomarkers in ensembles of parameterized boolean networks that have been subject to drug combination perturbations.

The `emba` R package aims to fill that gap and provide a first implementation of such a novel software.
Initially, it was designed as a complementary software tool, to help the analysis of the parameterized boolean model ensembles which were produced by modules from the DrugLogics NTNU software pipeline (see respective documentation [@dl-doc]).
Later, we generalized most of the functions in the package and modularized them to package-essential (that form the core of the `emba` package) and various general-purpose yet useful functions (that are now part of the dependency package `usefun` [@R-usefun]).

# Summary

The main functionality of the `emba` R package is to find *performance* and *synergy* biomarkers.
Performance biomarkers are nodes in the input boolean networks whose activity state and/or model parameterization affects the predictive performance of those models.
The prediction performance can be assessed via the number of true positive predictions or the Matthews correlation coefficient score which is more robust to class imbalances [@Chicco2020].
On the other hand, synergy biomarkers are nodes which provide hints for the mechanisms behind the complex process of synergy manifestation in drug combination datasets.

For more information, see our “Get started guide” and the reference manual in the package website [@emba-site].
Several analyses using the `emba` R package are available in a separate repository [@gma-git].
Future developments will include the implementation of a method for the identification of *topology* biomarkers, where we will be able to assess which interactions in the network are important for the manifestation of synergies in specific cell-contexts.

# Acknowledgements

This work was supported by ERACoSysMed grant *COLOSYS* (JZ, MK) and The NTNU Strategic Research Area *NTNU Health* (AF).

# References
