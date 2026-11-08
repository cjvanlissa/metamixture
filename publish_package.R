Sys.setenv(METAMIXTURE_RUN_EVERYTHING = "true")
Sys.setenv(METAMIXTURE_SHOW_RESULTS = "true")
f <- list.files(path = "vignettes", pattern = "Rmd$", full.names = TRUE)
for(thisf in f){
  rmarkdown::render(thisf, clean = FALSE)
}
Sys.setenv(METAMIXTURE_RUN_EVERYTHING = "false")
pkgdown::build_site()
Sys.setenv(METAMIXTURE_SHOW_RESULTS = "true")
devtools::check()
worcs::git_update("run everything in vignettes and update site")

