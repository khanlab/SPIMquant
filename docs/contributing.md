# Contributing to SPIMquant

Thank you for your interest in contributing to SPIMquant! This document provides guidelines for contributing to the project.

## Ways to Contribute

There are many ways to contribute to SPIMquant:

- **Report bugs**: Open issues for bugs you encounter
- **Suggest features**: Propose new features or improvements
- **Improve documentation**: Help us improve docs and examples
- **Submit code**: Fix bugs or implement new features
- **Share examples**: Contribute workflow examples and tutorials
- **Help others**: Answer questions in discussions

## Getting Started

### Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/SPIMquant.git
   cd SPIMquant
   ```
3. Install development environment:
   ```bash
   pixi install --environment dev
   ```
4. Create a branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

### Development Environment

The dev environment includes:
- Code formatters (black, isort)
- Snakemake formatter (snakefmt)
- Testing tools (pytest)
- Visualization tools (JupyterLab, napari)

## Code Standards

### Python Code Style

SPIMquant uses:
- **black** for code formatting
- **isort** for import sorting
- **snakefmt** for Snakemake files

Format your code before committing:

```bash
# Check formatting
pixi run quality_check

# Automatically fix formatting
pixi run quality_fix
```

### Coding Guidelines

1. **Follow PEP 8**: Python code should follow PEP 8 style guide
2. **Use type hints**: Add type hints to function signatures
3. **Write docstrings**: Document functions, classes, and modules
4. **Keep functions focused**: Each function should do one thing well
5. **Add tests**: Include tests for new functionality
6. **Comment complex code**: Explain non-obvious logic

### Snakemake Guidelines

1. **Modular rules**: Keep rules focused and reusable
2. **Use bids() function**: For BIDS-compliant file paths
3. **Document rules**: Add comments explaining rule purpose
4. **Parameter access**: Use `snakemake.params` in scripts
5. **Resource specification**: Declare memory/thread requirements

## Testing

### Running Tests

```bash
# Run all tests
pixi run pytest

# Run specific test file
pixi run pytest tests/test_workflow.py

# Run with verbose output
pixi run pytest -v
```

### Writing Tests

<!-- TODO: Add test writing guidelines when test infrastructure is established -->

1. Place tests in `tests/` directory
2. Follow existing test patterns
3. Test both success and failure cases
4. Use fixtures for common setup

### Dry Run Testing

Always test workflows with dry run:

```bash
pixi run spimquant /absolute/path/to/tests/bids_ds /tmp/output participant -n
```

## Making Changes

### Workflow for Changes

1. **Create an issue**: Describe the bug or feature
2. **Create a branch**: Use descriptive branch names
   - `feature/add-new-template-support`
   - `fix/registration-memory-leak`
   - `docs/improve-installation-guide`
3. **Make changes**: Follow code standards
4. **Test changes**: Run tests and dry runs
5. **Commit changes**: Use clear commit messages
6. **Push to fork**: Push your branch to GitHub
7. **Open pull request**: Submit PR for review

### Commit Messages

Write clear, descriptive commit messages:

```
Add support for custom templates

- Implement custom template loading
- Add configuration options
- Update documentation
- Add tests for custom templates

Closes #123
```

Format:
- First line: Brief summary (50 chars or less)
- Blank line
- Detailed explanation if needed
- Reference related issues

### Pull Request Guidelines

1. **Describe changes**: Explain what and why
2. **Link issues**: Reference related issues
3. **Include tests**: Add or update tests
4. **Update docs**: Document new features
5. **Follow template**: Use PR template if provided

## Documentation

### Improving Documentation

Documentation is in the `docs/` directory using MkDocs.

To build and serve locally:

```bash
# TODO: Add mkdocs commands once mkdocs is in dependencies
```

### Documentation Style

1. **Clear and concise**: Write for diverse audiences
2. **Use examples**: Include code examples
3. **Add TODOs**: Mark incomplete sections with `<!-- TODO: ... -->`
4. **Link related pages**: Help users navigate
5. **Test code examples**: Ensure examples work

### Adding New Documentation

1. Create new `.md` files in appropriate `docs/` subdirectory
2. Add to navigation in `mkdocs.yml`
3. Link from related pages
4. Follow existing structure and style

## Reporting Issues

### Bug Reports

When reporting bugs, include:

1. **Clear title**: Descriptive summary of the issue
2. **SPIMquant version**: Output of `pixi run spimquant --version`
3. **System information**: OS, memory, CPU details
4. **Steps to reproduce**: Minimal example to reproduce bug
5. **Expected behavior**: What should happen
6. **Actual behavior**: What actually happens
7. **Error messages**: Full error output and logs
8. **Additional context**: Anything else relevant

### Feature Requests

When suggesting features:

1. **Clear description**: What feature you'd like
2. **Use case**: Why it would be useful
3. **Proposed solution**: Ideas for implementation
4. **Alternatives**: Other options considered
5. **Additional context**: Examples from other tools

## Code Review Process

### What to Expect

1. **Review time**: May take a few days to weeks
2. **Feedback**: Reviewers may request changes
3. **Discussion**: Open dialogue about implementation
4. **Iteration**: May need multiple rounds of review
5. **Approval**: Maintainer approval required to merge

### Review Criteria

Reviews check for:

- Code quality and style
- Test coverage
- Documentation completeness
- Performance impact
- Breaking changes
- Security considerations

## Community Guidelines

### Code of Conduct

<!-- TODO: Add link to CODE_OF_CONDUCT.md if it exists -->

Be respectful and inclusive:

- Use welcoming language
- Respect differing viewpoints
- Accept constructive criticism gracefully
- Focus on what's best for the community
- Show empathy toward others

### Getting Help

- **Questions**: Use [GitHub Discussions](https://github.com/khanlab/SPIMquant/discussions)
- **Bugs**: Open [GitHub Issues](https://github.com/khanlab/SPIMquant/issues)
- **Chat**: <!-- TODO: Add chat/slack if available -->

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Recognition

Contributors are recognized in:
- GitHub contributors list
- Release notes
- Documentation acknowledgments

## Thank You!

Your contributions make SPIMquant better for everyone. Thank you for taking the time to contribute!