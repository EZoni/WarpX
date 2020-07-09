###### Make sure that the PR's description is up-to-date at every stage of the review process. If the PR is still a work in progress, convert the PR to a draft. If the PR is ready for review, mark the PR as ready for review. Please, refer to the GitHub Docs for more information about [changing the stage of a PR](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/changing-the-stage-of-a-pull-request). When your PR is ready for review, [request a review](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/requesting-a-pull-request-review) from both [WarpX admins](https://github.com/orgs/ECP-WarpX/teams/warpx-admins/members) and [WarpX contributors](https://github.com/orgs/ECP-WarpX/teams/warpx-contributors/members) and [assign the PR](https://docs.github.com/en/github/managing-your-work-on-github/assigning-issues-and-pull-requests-to-other-github-users) to at least one of the [WarpX admins](https://github.com/orgs/ECP-WarpX/teams/warpx-admins/members). Once you have addressed the reviewers' comments, please [re-request their review](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/requesting-a-pull-request-review) and continue to do so untill the PR is approved.

What is the goal of this PR?
---

Describe here the goal of your PR.

Description
---
###### You may [create permanent links](https://docs.github.com/en/github/managing-your-work-on-github/creating-a-permanent-link-to-a-code-snippet) of the most important code changes in order to guide reviewers through your new implementation.

Describe here the details of your PR.

How did I test this PR?
---
###### You may include relevant plots and tables as well as link to other PRs or issues where your implementation is merged and tested. You may also mention whether you added new automated tests to the CI test suite or modified existing tests, including the corresponding [checksum regression benchmarks](https://warpx.readthedocs.io/en/latest/developers/checksum.html).

Describe here how you tested your PR.

Checklists
--- 
###### Please check the appropriate items in the following lists to facilitate the review process.

Documentation:
* [ ] I added [Doxygen comments](https://www.doxygen.nl/manual/docblocks.html) (header files only)
* [ ] I updated and built the [online documentation](https://github.com/ECP-WarpX/WarpX/tree/development/Docs)

GPU portability:
* [ ] The GPU build is successful
* [ ] The GPU CI tests pass

Style guidelines:
* [ ] I used 4 spaces for indentation (no tabs)
* [ ] I kept the number of characters per line approximately < 100
* [ ] I added a space before and after all assignment operators `=`
* [ ] I added a space after function names (only function definitions)
* [ ] I did not include style changes that are not related to this PR
* [ ] I used the `CamelCase` naming convention for new files and classes
* [ ] I prefixed all the member variables of my new classes with `m_`
* [ ] I followed the guidelines described [here](https://github.com/ECP-WarpX/WarpX/blob/development/Docs/source/developers/repo_organization.rst) for `#include` directives
* [ ] I prefixed all AMReX types with `amrex::` and I included AMReX type literals with `using namespace amrex::literals`
