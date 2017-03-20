
#' R package for fitting temporal trends to abundence data aggregated over large regions when subregions have missing data
#' 
#' This package fits a log-linear trend models to regions aggregated over sites. The sites may contain missing surveys that 
#' are not temporally aligned with the missing data at other sites, making direct aggregation impossible. The functions within the package
#' model the indivdual sites with a semi-parametric (possibly, zero-inflated) model to interpolate missing data from which regional aggregations
#' can be made. By using Markov Chain Monte Carlo, on can sample from the posterior predictive distribution of the regional aggregations
#' Then calculate the log-linear trend over the time period of interest as a derived parameter. Using the posterior predictive distribution
#' allows incorporation of both parameter uncertainty as well as uncertainty due to sampling the local abundance processes.
#' 
#' \tabular{ll}{ Package: \tab agTrend\cr 
#' Type: \tab Package\cr 
#' Version: \tab 0.17.7\cr 
#' Date: \tab 2017-03-20\cr 
#' License: \tab Unlimited\cr 
#' LazyLoad: \tab
#' yes\cr }
#' 
#' @name agTrend-package
#' @aliases agTrend-package agTrend
#' @docType package
#' @author Devin S. Johnson
#' 
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#' 
#' 
NULL

#' Steller sea lion survey data (nonpup counts) collected by the NOAA National Marine Mammal Laboratory in the western Distinct 
#' Population Segment of the stock (wDPS).
#' 
#' 
#' @name wdpsNonpups
#' @docType data
#' @format
#' 
#' A data frame with 3408 observations on the following 12 variables.
#' 
#' \describe{ \item{site}{Site where survey was taken}
#' 
#' \item{Rookery}{1 if site is a rookery, 0 if else.}
#' 
#' \item{Region}{Region in which the site is located}
#' 
#' \item{RegNum}{Numeric representation of the 'Region' variable}
#' 
#' \item{RCA}{Rookery Cluster Area in which the site belongs}
#' 
#' \item{Lat}{Latitude of the site}
#' 
#' \item{Long}{Longitude of the site}
#' 
#' \item{trend70}{1 if the site belongs to the group of 70s trend sites}
#' 
#' \item{trend90}{1 if the site belongs to the group of 90s trend sites}
#' 
#' \item{CH}{Equals 1 if the site is within critical habitat boundries, 0 else.}
#' 
#' \item{year}{Year in which the survey was conducted}
#' 
#' \item{count}{The count of nonpups observed}
#' }
#' 
#' @references Need a tech report / memo reference...
#' 
#' @source Alaska Ecosystems Program, National Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' data(wdpsNonpups)
#' head(wdpsNonpups)
NULL

#' Data from a pilot study of collecting Steller sea lion counts from photos taken by hand at 
#' oblique angles vs. vertical medium-format photos.
#' 
#' 
#' @name photoCorrection
#' @docType data
#' @format
#' 
#' A data frame with 20 observations on the following 4 variables.
#' 
#' \describe{ \item{site}{Site where survey was taken}
#' 
#' \item{Region}{Region in which the site is located}
#' 
#' \item{X2000OBL}{Survey count taken from oblique photo source}
#' 
#' \item{X2000VERT}{Survey count taken from vertical medium-format photo source}
#' 
#' }
#' 
#' @source Alaska Ecosystems Program, National Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' data(photoCorrection)
#' head(photoCorrection)
NULL

#' Pup counts from surveys of the eastern Distict Populations Segment (eDPS) of Steller sea lions.
#'  
#' @name edpsPups
#' @docType data
#' @format
#' 
#' A data frame with 35 observations on the following 3 variables.
#' 
#' \describe{ 
#' \item{region}{One of 4 eDPS subregions: southeast Alaska (SE AK), British Columbia (BC), Oregon (OR), and California (CA)}
#' 
#' \item{year}{Year the survey was conducted}
#' 
#' \item{count}{Aggregated count of pups of sites within each subregion.}
#' }
#' 
#' @source Alaska Ecosystems Program, National Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' data(edpsPups)
#' head(edpsPups)
NULL

#' Pup counts from surveys of the western Distict Populations Segment (wDPS) of Steller sea lions.
#'  
#' @name wdpsPups
#' @docType data
#' @format
#' 
#' A data frame with 35 observations on the following 3 variables.
#' 
#' \describe{ 
#' \item{site}{Surveyed sites in the wDPS}
#' 
#' \item{Region}{Region of the surveyed site}
#' 
#' \item{RCA}{Rookery cluster area of the surveyed site}
#' 
#' \item{year}{Survey year}
#' 
#' \item{count}{Survey count of pups at each site}
#' }
#' 
#' @source Alaska Ecosystems Program, National Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' data(wdpsPups)
#' head(wdpsPups)
NULL



.onAttach <- function(library, pkgname)
  {
    info <-utils::packageDescription(pkgname)
    package <- info$Package
    version <- info$Version
    date <- info$Date
    packageStartupMessage(
      paste("\n\n",paste(package, version, paste("(",date, ")", sep=""), "\n\n"), 
            "Type 'demo(package='agTrend')' to see a list of demos for this package.\n\n",
            "The raw code for the demos can be found by typing:\n",
            "\t system.file('demo', package='agTrend')\n\n",
            "To access the help files type:\n",
            "\t help('agTrend-package')\n\n"
            )
    )
  }

