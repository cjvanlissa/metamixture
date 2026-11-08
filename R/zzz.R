.onLoad <- function(libname, pkgname) {

  tidySEM::register_tidySEM_method(
    generic = "random_starts_internal",
    type = "remix",
    fun = random_starts_internal_remix
  )
  tidySEM::register_tidySEM_method(
    generic = "extract_postprob",
    type = "remix",
    fun = extract_postprob_remix
  )
  tidySEM::register_tidySEM_method(
    generic = "calc_fitindices",
    type = "remix",
    fun = calc_fitindices_mixture
  )
  tidySEM::register_tidySEM_method(
    generic = "two_step_internal",
    type = "remix",
    fun = two_step_internal_remix
  )
  tidySEM::register_tidySEM_method(
    generic = "num_classes",
    type = "remix",
    fun = num_classes_remix
  )
  tidySEM::register_tidySEM_method(
    generic = "predict_class",
    type = "remix",
    fun = predict_class_remix
  )
  tidySEM::register_tidySEM_method(
    generic = "auxiliary_internal",
    type = "remix",
    fun = auxiliary_internal_remix
  )

}
