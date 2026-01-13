# Contributing to mvndist

Thank you for considering contributing to mvndist! 

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:

- A clear description of the problem
- A minimal reproducible example
- Your R version and operating system
- Expected vs. actual behavior

### Suggesting Enhancements

Feature requests are welcome! Please open an issue describing:

- The use case for the feature
- How it would work
- Any relevant references or implementations

### Pull Requests

1. Fork the repository
2. Create a new branch for your feature (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Update documentation (roxygen2 comments)
6. Run `devtools::check()` to ensure everything passes
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

## Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable and function names
- Add roxygen2 documentation for all exported functions
- Include examples in documentation

## Testing

- Add unit tests for new functions in `tests/testthat/`
- Ensure all tests pass with `devtools::test()`
- Check code coverage with `covr::package_coverage()`

## Documentation

- Update roxygen2 comments for any function changes
- Run `devtools::document()` to regenerate `.Rd` files
- Update README.md if adding major features
- Add entries to NEWS.md

## Release communication (R-bloggers)

- For any CRAN release or major feature, draft a short blog post summarizing what changed.
- Include 3-5 bullet highlights, installation instructions, and one minimal runnable example.
- Publish the post on a blog feed registered with R-bloggers; if the feed is not registered yet, submit it once via the R-bloggers site.
- Link back to the package repo, NEWS.md, and (if available) the pkgdown site; note the package version and release date.
- Update or append the post if follow-up fixes significantly change usage or behavior.

## Questions?

Feel free to open an issue for questions about contributing!
