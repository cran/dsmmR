---
title: NEWS
editor_options: 
  markdown: 
    wrap: 72
---

# dsmmR 1.0.5

## Tests

- Tests are now implemented for the functions of the package. Namely:
    - `fit_dsmm()`
    - `create_sequence()`
    - `get_kernel()`
    - `parametric_dsmm()`
    - `nonparametric_dsmm()`
  (@Bisaloo, 2)
    
## Code of Conduct

- The Contributor Covenant Code of Conduct is now added with respect to the JOSS 
  publication. The `CODE_OF_CONDUCT.md` file can be found in the `.github` folder.

## Documentation

- `get_kernel()` has small changes in the description, regarding the
  reduction of dimensions when selecting a specific argument of u, v, l or t.

## Workflows

- The Github workflow `R-hub v2` is added in `.github/workflows`. This ensures
  that after every github update of the online `dsmmR` repository, the CRAN
  checks are being made for multiple platforms (linux, macos, windows). This can
  be run manually through the `command rhub::rhub_check(branch = 'master')`.
  
- Another Github workflow `codecov` was added in `.github/workflows`. This 
  ensures that the `dsmmR` package is mostly covered by the automated tests in 
  place.

## Errors

- Fixed an unexpected check that was failing in the functions `valid_p_dist()`,
  `valid_fdist_nonparametric()`, `valid_fdist_parametric()` and added a cases for
  the f non-drifting case in `get_fdist_parametric()`. This was causing some
  errors to appear when trying to print an estimated fitted model when f was not
  drifting (Model 2).

## Reference DOIs

- For each reference cited in `paper.bib`, DOIs were added in the correct
  formatting. This was a small fix done for the JOSS publication.

## Error handling

- Increased the clarity of the string formatting for the different warning
  messages the user may encounter during his use of the package.

# dsmmR 1.0.4

## DESCRIPTION

- `fit_dsmm()` gains a new attribute `multi_estimation`, which enables the 
  estimation of a drifting semi-Markov model using multiple sequences. There 
  are two possible options: `avg_model` and `count_sum`.
     - `avg_model` averages the `q_i` received from multiple sequences 
     - `count_sum` adds the counts of the states for each sequence
       (of equal size) and then computes the `q_i`. 
       
- `simulate.dsmm()` gains a new attribute, `max_seq_length`, which is renamed
  from old attribute `seq_length` for clarity.
  

# dsmmR 1.0.3

## DESCRIPTION

- Updated the `Depends` section. Now we impose the requirement for R >= 3.5.0,
  in order to make proper use of the `isTRUE()` and `isFALSE()` functions in 
  the `is_logical()` function defined in `utils.R`. These functions will remain
  for their clarity.
  

## Documentation

- `fit_dsmm()` and `simulate.dsmm()` functions now better explain the difference
  between a sequence of states and the embedded Markov chain.
  Notably, it was specified that sequences are input in `fit_dsmm()`, specified
  with the argument `sequence`.
  In `simulate.dsmm()`, such a sequence is the resulting output. 
  The embedded Markov chain can be seen as part of the output from `fit_dsmm()`,
  named `emc`. Furthermore, the function `base::rle()` is mentioned for clarity.


# dsmmR 1.0.2

## Documentation

- Updated the `README` and `DESCRIPTION` files with an acknowledgement section.


# dsmmR 1.0.1

## Minor Improvements

-   Added a `NEWS.md` file to track changes to the package.

-   Now the `fit_dsmm()` function has a default value for the `states` attribute,
    being the sorted unique values of the `sequence` character vector attribute.
    

## Bug fixes

-   Fixed a case where `simulate.dsmm()` sometimes did not function as expected
    when `nsim = 1`.
    -   Now it is possible to specify `nsim = 0`, so that the simulated
        sequence will only include the initial state and its
        corresponding sojourn time, e.g. "a", "a", "a".
        
        By giving `nsim = 1` , a single simulation will be made from the
        drifting semi-Markov kernel, returning for example "a", "a",
        "a", "c".

## Documentation

-   Updated the documentation for `simulate.dsmm()`, with accordance to
    the changes made.

-   Updated the `README` file.

    -   Added high-level documentation of the package.
    
    -   Added installation instructions with access to
        the development version of the package through github.

-   Updated the documentation for `dsmmR-package`.

    -   Added a "Community Guidelines" section, so that users can report
        errors or mistakes and contribute directly to the software
        through the newly-established open-source github page at
        <https://github.com/Mavrogiannis-Ioannis/dsmmR>.

    -   Added a "Notes" section, specifying that automated tests are in
        place in order to aid the user with any false input made and,
        furthermore, to ensure that the functions used return the
        expected output.
