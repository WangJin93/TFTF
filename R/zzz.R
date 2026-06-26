#' @importFrom utils packageDescription

.onAttach <- function(libname, pkgname) {
  # 获取包版本信息
  pkg_info <- packageDescription("GCAS")
  version <- pkg_info$Version
  
  # 构建欢迎消息
  welcome_msg <- c(
    paste0("GCAS v", version, " loaded successfully!"),
    "",
    "Citation:",
    "  Wang, J.; Wei, M.; Zhang, J.; Song, X.; Hu, Y.; Qin, L.; Liang, T.; Zhu, X.; Li, J.",
    "  GCAS: An Integrated R Package and Shiny App for Comprehensive Cancer Data Analysis.",
    "  Biomolecules 2026, 16, 823.",
    "",
    "For more information, visit: https://wangjin93.github.io/gcas.html",
    "Run GCAS::GCAS_app() to launch the Shiny app."
  )
  
  # 使用 packageStartupMessage 输出消息
  packageStartupMessage(paste(welcome_msg, collapse = "\n"))
}
