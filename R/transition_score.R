
#' Computes the local transition score
#'
#' This function computes the local transition score, i.e. a score between 0 and 1 (non-linearly) correlated with the distance to a transition from or to a CST-X sample.
#' @param transition_type a character, either \code{"IN"} or \code{"OUT"} which indicates if the score has to be computed for transitions TO ("IN") CST4 or FROM ("OUT") CST4.
#' @param CST_seq specifies the CST sequence for which the transition score should be computed.
#'      This can be specified as a vector of
#'      (1) integers in 1:5 representing the CSTs or
#'      (2) characters sampled from \code{c("CST1", "CST2", "CST3", "CST4", "CST5")} or
#'      (3) logicals in which \code{TRUE} represent samples with the CST of interest.
#' @param time_seq specifies the time-points at which the samples were taken. Should be specified as a numerical vector.
#' @param CST_X (optional) (default = 4 or "CST4") specifies the CST of interest. If \code{CST_seq} is a vector of logicals, the value of this argument is irrelevant.
#' @param timescale (optional) a strictly positive number (double) which specifies the time-scale of the transition score. Must be expressed in the same units as the time_seq. Default value is 4.
#' @param kernel_fun a function which specifies the kernel that should be use to transform the distance to the transition to the transition score.
#'      The default is a sigmoid kernel.
#'      The specified function must have \code{'time'} and \code{'scale'} as arguments.
#'
#' @return a vector of scores between 0 and 1 and of NAs for time-points where the score could not be computed.
#'
#' @export
#' @importFrom magrittr %>%
#' @examples
#'

compute_local_transition_score = function(
  data = data.frame(),
  transition_type = c("IN","OUT"),
  CST_X = 4,
  timescale = 4,
  kernel_fun = shifted_scaled_sigmoid
){

  # CHECKS
  if(nrow(data) == 0) return(data)
  data_required_cols = c("subject_id","t","CST")
  if(!all(data_required_cols %in% colnames(data))) stop(paste0("'data' must have the following columns: ",paste0(data_required_cols, collapse = ", ")))
  if(!is.numeric(data$t)) stop("column 't' of 'data' must be numeric.")
  if(!is.integer(data$CST)) stop("column 'CST' of 'data' must be composed of integer.")


  transition_type = transition_type[1]
  if(! transition_type %in% c("IN","OUT")) stop("'transition_type' must be 'IN' or 'OUT'.")
  if(timescale < 0 ) stop("'timescale' must be a positive number.")

  if((CST_X %% 1)!= 0) stop("'CST_X' must be an integer.")
  CST_X = CST_X[1]

  # order the sequence
  data = data %>% dplyr::arrange(subject_id, t)

  # define x:
  data = data %>% dplyr::mutate(x = (CST == CST_X))
  if(transition_type == "OUT") data = data %>% dplyr::mutate(x = !x)

  # identify "off" transitions
  data = data %>%
    dplyr::group_by(subject_id) %>%
    dplyr::mutate(transition = (!x & dplyr::lag(x)) %>% tidyr::replace_na(FALSE),
                  transition_id = cumsum(transition))

  # identify time of next transition and time to transition
  data = data %>%
    dplyr::group_by(subject_id, transition_id) %>%
    dplyr::mutate(next_transition_t = suppressWarnings(min(t[x])),
                  time_to_transition = (t - next_transition_t) %>% pmin(0),
                  time_to_transition = ifelse(time_to_transition == 0, NA, time_to_transition)) %>%
    dplyr::ungroup()

  # set time-points that are within timescale of end of time-series to NAs
  data = data %>%
    dplyr::group_by(subject_id) %>%
    dplyr::mutate(last_t = max(t),
                  close_to_end = (last_t - t) < timescale,
                  time_to_transition = ifelse(close_to_end & is.infinite(time_to_transition), NA, time_to_transition)) %>%
    dplyr::ungroup()



  # compute_score
  data = data %>%
    dplyr::mutate(score = kernel_fun(t = time_to_transition, timescale = timescale, tx = -0.5))

  # clean the data.frame
  data = data %>% dplyr::select(subject_id, t, score)
  # return results
  data
}


#' #' Visualizes local transition score
#' #'
#' #' This function visualizes the local transition score.
#' #' @seealso compute_local_transition_score
#' #' @param df a \code{data.frame} with the following columns: subject_it (num or character or factor), t (time-point, must be numeric), CST (integer in 1:5), Transition_score_IN, Transition_score_OUT
#' #'
#' #' @return a \code{ggplot}.
#' #'
#' #'
#' #' @import ggplot2
#'
#' plot_transition_scores = function(df){
#'
#'   g = ggplot(df, aes(x = Timepoint_week))+
#'     geom_line(aes(y = Transition_score_IN), col = "tomato")+
#'     geom_point(aes(y = Transition_score_IN), col = "tomato")+
#'     geom_line(aes(y = Transition_score_OUT), col = "steelblue")+
#'     geom_point(aes(y = Transition_score_OUT), col = "steelblue")+
#'     geom_point(aes(y = 1.5, col = factor(CST_num)),size = 3, alpha = 0.6)+
#'     scale_color_manual(values = c("steelblue1","steelblue2","steelblue3","red","steelblue4","steelblue"))+
#'     facet_wrap(SID ~ . )+
#'     theme_set(theme_minimal())
#'   g
#'
#' }




sigmoid = function(t, slope, t0){
  x = 1 / (1 + exp(-slope*(t - t0)))
  x
}


shifted_scaled_sigmoid = function(t, timescale, tx){
  margin = 0.025
  x = sigmoid(t, slope =  2 * log((1-margin)/margin) / timescale, t0 = tx - timescale/2)
}

