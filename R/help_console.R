#' Print help file into console
#'
#' From https://stackoverflow.com/questions/12038160/how-to-not-run-an-example-using-roxygen2
#'   based on code by Noam Ross
#'   http://www.noamross.net/archives/2013-06-18-helpconsoleexample/
#'   Stéphane Laurent
#'   https://stackoverflow.com/questions/60468080/
#'   print-an-r-help-file-vignette-as-output-into-an-r-html-notebook
#'   and Michael Sumner (mdsumner)
#'   https://stackoverflow.com/questions/7495685/
#'   how-to-access-the-help-documentation-rd-source-files-in-r
#'
#' @param topic the command for which help is required
#' @param package the package name with the required topic
#' @param format output format
#' @param before place code before the output e.g. `<blockquote>`
#' @param after place code after the output e.g. `</blockquote>`
#'
#' @return formatted function help
#'
#' @export
help_console <- function(topic,
                         package = "replicationOrigins",
                         format = c("text", "html", "latex", "Rd"),
                         before = NULL, after = NULL) {

  prepare_Rd <- utils::getFromNamespace("prepare_Rd", "tools")
  format <- match.arg(format)
  if (!is.character(topic)) topic <- deparse(substitute(topic))
  db <- tools::Rd_db(package)
  helpfile <- db[paste0(topic, ".Rd")][[1]]

  hs <- utils::capture.output(
    switch(
      format,
      text = tools::Rd2txt(helpfile),
      html = tools::Rd2HTML(
        helpfile,
        package = "",
        stages = c("install", "render")
      ),
      latex = tools::Rd2latex(helpfile),
      Rd = prepare_Rd(helpfile)
    )
  )
  if (format == "html") {
    i <- grep("<body>", hs)
    j <- grep("</body>", hs)
    hs <- hs[(i+1):(j-1)]
  }
  hs <- stringr::str_replace_all(hs, "h2>", "h3>")
  hs <- stringr::str_replace_all(hs, "h3>", "h4>")
  hs <- c(before, hs, after)
  hs <- cat(hs, sep = "\n")
  invisible(hs)
}


