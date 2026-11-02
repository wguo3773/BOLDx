pkg <- normalizePath(".")
roxygen2::roxygenise(pkg)
for (nm in c("bold_config","simulate_bold","simulate_scenarios","summarise_scenarios",
             "beta_mixture","recommend_dose","plan_intervention")) {
  rd <- tools::parse_Rd(file.path(pkg,"man",paste0(nm,".Rd")))
  tools::Rd2txt(rd, out=file.path(pkg,"docs",paste0(nm,".txt")),
    options=list(underline_titles=FALSE))
}
