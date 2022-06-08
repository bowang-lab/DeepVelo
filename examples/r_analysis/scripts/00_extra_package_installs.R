require(devtools)

# Install R packages not available in the Conda repository
install_version(
    "ActivePathways",
    version = "1.1.0",
    repos = "http://cran.us.r-project.org"
)
install.packages(
    "ggvenn",
    repos = "http://cran.us.r-project.org"
)