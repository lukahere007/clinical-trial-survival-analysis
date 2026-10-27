# Load the project's helper functions (tests run from tests/testthat).
for (f in list.files(file.path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
