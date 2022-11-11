
#' The Bone Marrow Transplant Data
#' 
#' Bone marrow transplant data with 408 rows and 5 columns.
#' 
#' 
#' @format The data has 408 rows and 5 columns. \describe{ \item{cause}{a
#' numeric vector code.  Survival status. 1: dead from treatment related
#' causes, 2: relapse , 0: censored.} \item{time}{ a numeric vector. Survival
#' time.  } \item{platelet}{a numeric vector code. Plalelet 1: more than 100 x
#' \eqn{10^9} per L, 0: less.} \item{tcell}{a numeric vector. T-cell depleted
#' BMT 1:yes, 0:no.} \item{age}{a numeric vector code. Age of patient, scaled
#' and centered ((age-35)/15).} }
#' @references NN
#' @name bmt
#' @docType data
#' @source Simulated data
#' @keywords package 
#' @examples
#' 
#' data(bmt)
#' names(bmt)
#' 
NULL


#' The Diabetic Retinopathy Data
#' 
#' The data was colleceted to test a laser treatment for delaying blindness in
#' patients with dibetic retinopathy. The subset of 197 patiens given in Huster
#' et al. (1989) is used.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{id}{a numeric vector. Patient code.} \item{agedx}{a numeric vector.
#' Age of patient at diagnosis.} \item{time}{a numeric vector. Survival time:
#' time to blindness or censoring.} \item{status}{ a numeric vector code.
#' Survival status. 1: blindness, 0: censored.} \item{trteye}{a numeric vector
#' code. Random eye selected for treatment. 1: left eye 2: right eye.}
#' \item{treat}{a numeric vector. 1: treatment 0: untreated.} \item{adult}{a
#' numeric vector code. 1: younger than 20, 2: older than 20.} }
#' @source Huster W.J. and Brookmeyer, R. and Self. S. (1989) MOdelling paired
#' survival data with covariates, Biometrics 45, 145-56.
#' @name  diabetes
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(diabetes)
#' names(diabetes)
#' 
NULL


#' The Melanoma Survival Data
#' 
#' The melanoma data frame has 205 rows and 7 columns.  It contains data
#' relating to survival of patients after operation for malignant melanoma
#' collected at Odense University Hospital by K.T.  Drzewiecki.
#' 
#' 
#' @format This data frame contains the following columns: \describe{
#' \item{no}{ a numeric vector. Patient code. } \item{status}{ a numeric vector
#' code. Survival status. 1: dead from melanoma, 2: alive, 3: dead from other
#' cause. } \item{days}{ a numeric vector. Survival time. } \item{ulc}{ a
#' numeric vector code. Ulceration, 1: present, 0: absent. } \item{thick}{ a
#' numeric vector. Tumour thickness (1/100 mm). } \item{sex}{ a numeric vector
#' code. 0: female, 1: male. } }
#' @source Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993),
#' \emph{Statistical Models Based on Counting Processes}, Springer-Verlag.
#' 
#' Drzewiecki, K.T., Ladefoged, C., and Christensen, H.E. (1980), Biopsy and
#' prognosis for cutaneous malignant melanoma in clinical stage I. Scand. J.
#' Plast. Reconstru. Surg. 14, 141-144.
#' @name  melanoma 
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(melanoma)
#' names(melanoma)
#' 
NULL





#' The TRACE study group of myocardial infarction
#' 
#' The TRACE data frame contains 1877 patients and is a subset of a data set
#' consisting of approximately 6000 patients.  It contains data relating
#' survival of patients after myocardial infarction to various risk factors.
#' 
#' sTRACE is a subsample consisting of 300 patients.
#' 
#' tTRACE is a subsample consisting of 1000 patients.
#' 
#' 
#' @aliases TRACE sTRACE tTRACE
#' @format This data frame contains the following columns: \describe{
#' \item{id}{a numeric vector. Patient code. } \item{status}{ a numeric vector
#' code. Survival status. 9: dead from myocardial infarction, 0: alive, 7: dead
#' from other causes.  } \item{time}{ a numeric vector. Survival time in years.
#' } \item{chf}{ a numeric vector code. Clinical heart pump failure, 1:
#' present, 0: absent. } \item{diabetes}{ a numeric vector code. Diabetes, 1:
#' present, 0: absent. } \item{vf}{ a numeric vector code. Ventricular
#' fibrillation, 1: present, 0: absent. } \item{wmi}{ a numeric vector.
#' Measure of heart pumping effect based on ultrasound measurements where 2 is
#' normal and 0 is worst. } \item{sex}{ a numeric vector code. 1: female, 0:
#' male. } \item{age}{ a numeric vector code. Age of patient. } }
#' @source The TRACE study group.
#' 
#' Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L., F. E. Nielsen,
#' Melchior, T., Joen, T. and P. K. Andersen (1997), Does in-hospital
#' ventricular fibrillation affect prognosis after myocardial infarction?,
#' European Heart Journal 18, 919--924.
#' @name  TRACE
#' @docType data
#' @keywords package 
#' @examples
#' 
#' data(TRACE)
#' names(TRACE)
#' 
NULL



