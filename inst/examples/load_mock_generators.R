# Load mock data generators
# This script loads the mock data generator functions for speciateIT

# Find the full path to the mock_generators.R file
mock_generators_path <- system.file("examples", "mock_generators.R", package = "microFGT")

# If the file doesn't exist in the package, look in the development directory
if (mock_generators_path == "") {
  # Try the development directory
  mock_generators_path <- "R/data/mock_generators.R"
  if (!file.exists(mock_generators_path)) {
    stop("Could not find mock_generators.R file. Make sure you're running this from the package root directory.")
  }
}

# Source the file
source(mock_generators_path)

# Let the user know it worked
message("Mock data generator functions loaded successfully.")
message("Available functions:")
message(" - generate_mock_speciateit()")
message(" - write_mock_speciateit()")
message(" - convert_speciateit_to_fgt()")