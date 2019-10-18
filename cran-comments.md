## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

- Please omit the redundant "A set of functions that help ... " and start with sth like "Analyses and visualizes ...".

I changed the `Description` accordingly.

- Please make sure that you do not change the user's options, par or working directory. If you really have to do so, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited:
```
oldpar <- par(mfrow=c(2,2))             # code line i
on.exit(par(oldpar))                    # code line i + 1
```

I don't change the user's options or working directory anywhere in the code 
and the 3 uses I had of `par(mar)`, I fixed them as proposed.

- Please also explain this in the description as most ppl wouldn't know. (referring to the use of the word 'DrugLogics' in the description)

I did so.

- Please add references describing the methods in your package to the description field of your DESCRIPTION file in the form:

authors (year) <doi:...>  
authors (year) <arXiv:...>  
authors (year, ISBN:...)  
or only if none those are available:  <https:...>  
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.

Currently the DrugLogics software pipeline paper is a draft and not citable (I could add a paper/doi URL later I guess?). 

The methods I used are explained well enough I believe in each respective function (I have added sufficient 'Description', 'Arguments' and 'Details' fields). I have an analysis that is used to showcase the package's functions ([see it here](https://bblodfon.github.io/gitsbe-model-analysis/atopo/cell-lines-2500/)), whose link I had included in the `README.md` and that's why I didn't include any other URLs in the DESCRIPTION (unless you want me to).
