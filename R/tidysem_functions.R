# #' @importFrom tidySEM table_fit mx_diagnostics
# #' @export
# tidySEM::table_fit
#
# #' @importFrom tidySEM table_results
# #' @export
# tidySEM::table_results

#' @importFrom utils getFromNamespace
is_remix <- utils::getFromNamespace("is_remix", "tidySEM")
num_classes <- utils::getFromNamespace("num_classes", "tidySEM")
entropy <- utils::getFromNamespace("entropy", "tidySEM")
label_auxiliary_parameters <- utils::getFromNamespace("label_auxiliary_parameters", "tidySEM")
calc_fitindices_mixture <- utils::getFromNamespace("calc_fitindices_mixture", "tidySEM")
extract_postprob <- utils::getFromNamespace("extract_postprob", "tidySEM")
label_auxiliary_parameters <- utils::getFromNamespace("label_auxiliary_parameters", "tidySEM")
mx_starts_auxiliary <- utils::getFromNamespace("mx_starts_auxiliary", "tidySEM")
auxiliary_internal <- utils::getFromNamespace("auxiliary_internal", "tidySEM")
