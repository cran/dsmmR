---
title: NEWS
editor_options: 
  markdown: 
    wrap: 72
---

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

    -   Added high-level documentation of the package 
    
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
