# Contribution guide

## How to start?

You can contribute in multiple ways:
- Reporting an issue
- Requiring a new feature
- Testing the code and the software in different PC configurations and OS
- Improving the code
- Checking the documentation content
- Generating a new or improving the current BFE methods
- and much more...

For issue reporting we have created a template, which contains almost everything needed to identify the problem.

## Step-by-step guide

### Before reporting an issue, requesting a feature, or asking a question
Please ensure that you have read the following docs:
- [documentation and FAQ](docs/Q&A/README.md)
- [minimal examples section](examples/README.md)
- [previous reported issues](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues)
- [previous discussions](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/discussions)

### New feature

1. Open an issue with your feature description;
2. We shall discuss the design and its implementation details;
3. Once we agree that the plan looks good, go ahead and implement it.


### Bugfix

1. Go to [GitHub issues](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues);
2. Pick an issue and comment on the task that you want to work on this feature;
3. If you need more context on a specific issue, please ask, and we will discuss the details.

Once you finish implementing a feature or bugfix, please send a Pull Request.

If you are not familiar with creating a Pull Request, here are some guides:
- http://stackoverflow.com/questions/14680711/how-to-do-a-github-pull-request
- https://help.github.com/articles/creating-a-pull-request/


## Contribution best practices

- Keep scientific/model changes separate from documentation and formatting changes.
- Add or update focused tests for changed behavior and document any limitations.
- Keep canonical example READMEs under `examples/` authoritative; run
  `python scripts/sync_example_docs.py --check` after editing them.
- Run the relevant focused tests and `mkdocs build --strict` before opening a pull request.


## GitHub CI

The continuous-integration workflow runs the project's test and documentation checks.
Reproduce the documentation gate locally with `mkdocs build --strict` after installing
`docs/requirements.txt`.

### Documentation

The gmx_MMPBSA's documentation is based in Mkdocs-Material. We additionally use external plugins and modified code to optimize the content.

If you have some issues with building docs - please make sure that you installed the required pip packages.

### Tests

Run the full unit-test discovery when practical:

```bash
python -m unittest discover -s tests
```

For a focused change, run the affected test module and record the exact command and
result in the pull request.

#### Adding new tests

Add regression coverage for user-visible behavior, parser contracts, documentation
examples, or manifest changes. Keep fixtures small and avoid silently changing
scientific settings to make a test pass.

### Integrations

If a contribution adds dependencies or an external workflow, update the relevant
installation/environment documentation and add a focused validation path. Keep the
supported dependency ranges aligned with `setup.py`, `docs/env.yml`, and the release
documentation.
