% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Load.r.functions.R
\name{get_mol}
\alias{get_mol}
\title{Get metabolites whose annotated type belongs to molecule.}
\usage{
get_mol(colon_sub)
}
\arguments{
\item{colon_sub}{A subset of 'colon_cancer_desi.clustered_mass.table.with.anno.csv'}
}
\value{
Four matrices with different levels of screening.
        mol_T: Molecule annotations are not empty.
        mol_T2: The names of annotations that are not all 'NONE'.
        mol_T2uni: The annotations is non-empty, non-NONE, and unique.
}
\description{
Get metabolites whose annotated type belongs to molecule.
}
\examples{
" "

}
