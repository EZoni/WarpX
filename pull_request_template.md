###### Make sure that your PR's description is up-to-date at every stage of the review process: if your PR is still a work in progress, convert it to a draft; if your PR is ready for review, mark it as ready for review (learn [here](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/changing-the-stage-of-a-pull-request) how to change the stage of your PR). When your PR is ready for review, [request a review](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/requesting-a-pull-request-review) from both [WarpX admins](https://github.com/orgs/ECP-WarpX/teams/warpx-admins/members) and [WarpX contributors](https://github.com/orgs/ECP-WarpX/teams/warpx-contributors/members) and [assign the PR](https://docs.github.com/en/github/managing-your-work-on-github/assigning-issues-and-pull-requests-to-other-github-users) to at least one of the [WarpX admins](https://github.com/orgs/ECP-WarpX/teams/warpx-admins/members). Once you have addressed all reviewers' comments, [request a new review](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/requesting-a-pull-request-review) and continue to do so until your PR is approved.

What is the goal of your PR?
---

Describe here the goal of your PR.

Description
---
###### In order to guide reviewers through your implementation, you may [create permanent links](https://docs.github.com/en/github/managing-your-work-on-github/creating-a-permanent-link-to-a-code-snippet) to a specific line or range of lines of code highlighting the most important changes introduced in your PR.

Describe here the details of your PR.

How did you test your PR?
---
###### You may include relevant plots and tables as well as link to other PRs or issues where your implementation is merged and tested. You may also mention whether you added new automated tests to the CI test suite or modified existing tests, including the corresponding [checksum regression benchmarks](https://warpx.readthedocs.io/en/latest/developers/checksum.html).

Describe here how you tested your PR.

Checklists
--- 
###### You may check the appropriate items in the following lists in order to facilitate the review process.

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
