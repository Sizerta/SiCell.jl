# Contributing to SiCell.jl

Thank you for your interest in SiCell.jl! Contributions of all kinds are welcome: bug reports, feature requests, documentation improvements, tutorials, and code.

## Reporting bugs

Please open an issue on the [GitHub issue tracker](https://github.com/Sizerta/SiCell.jl/issues) and include:

- a short description of the problem and what you expected to happen,
- a minimal reproducible example (ideally a small dataset or a snippet that runs on simulated data),
- your Julia version (`versioninfo()`) and SiCell.jl version (`] status SiCell`),
- the full error message and stack trace, if any.

## Requesting features

Open an issue describing the analysis you want to do, why existing functionality doesn't cover it, and, if possible, a reference to the method (paper, other package) you have in mind. Discussing a feature in an issue before writing code helps make sure it fits the package's design.

## Contributing code

1. Fork the repository and create a branch from `main` with a descriptive name (e.g. `fix-normalization-sparse`, `add-trajectory-plot`).
2. Set up a development environment:
   ```julia
   ] dev path/to/your/SiCell.jl
   ] test SiCell
   ```
3. Make your changes. Please:
   - follow the existing code style and naming conventions,
   - add docstrings to new exported functions,
   - add tests under `test/` for new functionality or bug fixes.
4. Run the full test suite locally (`] test SiCell`) and make sure it passes.
5. If your change affects user-facing behaviour, update the documentation and add an entry to `CHANGELOG.md`.
6. Open a pull request against `main`, describing what the change does and linking any related issue. Continuous integration will run the tests automatically.

Small, focused pull requests are easier to review than large ones. If you're planning a substantial change, please open an issue first.

## Improving documentation

Documentation fixes — typos, unclear explanations, new examples or tutorials — are very valuable and are a great first contribution. The same pull-request process applies.

## Support and maintenance

SiCell.jl is maintained by Sizerta (ORCID: [0000-0002-8602-8092](https://orcid.org/0000-0002-8602-8092)).

- **Questions and usage help:** open an issue with the `question` label.
- **Response times:** the maintainer aims to respond to issues and pull requests within two weeks.
- **Releases:** new versions are tagged on GitHub, registered in the Julia General registry, and archived on Zenodo. Changes are listed in `CHANGELOG.md`.
- **Decisions:** design and release decisions are currently made by the maintainer, taking discussion in issues and pull requests into account. Regular contributors may be invited to become co-maintainers.

## Code of conduct

Please be respectful and constructive in all interactions. We follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/) code of conduct.

## License

By contributing, you agree that your contributions will be licensed under the same license as SiCell.jl (see `LICENSE`).
