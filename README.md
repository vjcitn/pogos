# Overview

[PharmacoDB](https://pharmacodb.ca) unites multiple drug sensitivity databases.
A REST API is available and interfaces to some basic endpoints are
provided in this package.  We also provide some basic
support for working with the Cell Line Ontology and ChEBI.

A basic purpose of this package is to support exploration of
the completeness and utility 
of available ontologies for cell lines and drugs.
A relevant study of the possibility of "drug set
enrichment analysis" has been published
by [Napolitano and colleagues](https://www.ncbi.nlm.nih.gov/pubmed/26415724).

We note that detailed RESTful interrogation of many ontologies
is supported by the `r Biocpkg("rols")` package.  In the long
run we would like to make use of the EBI Ontology Lookup System
as a central resource of this package, but to develop concepts
at these early stages, static images of ontologies are employed.

# Installation

pogos can be installed using Bioconductor.  Be sure
BiocManager has been installed from CRAN.  Then:
```
BiocManager::install("pogos")
```
will take care of all dependencies.

Use the "Get started" tab above to learn more.
