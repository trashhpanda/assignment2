# Assignment 2: Bootstrap

## Setup

So that I can save myself the 2h of figuring it out for the next assignment.

1. Add the `.Rnw` template to the project
2. Run `renv::init()` then choose `2` to discard the default lockfile and initialize the project with installing the needed packages.
3. After that, `renv::status()` should output "No issues found -- the project is in a consistent state."
4. Try to compile the report template with `knitr`, see section "Compiling the report template".
5. Then just to keep it for reference, copy the template and rename the new file as `report.pdf`. That's the one you'll be editing.
6. Add a `source.R` to have a safe space for coding that doesn't require you to recompile the entire pdf each time you want to see if your code works. 


## Compiling the report template

After the setup, `knitr` should already be there. What's missing is `tinytex` so:

```
install.packages("tinytex")
tinytex::install_tinytex()
```

Then just `knitr::knit2pdf("template.Rnw")` should produce a pdf.

Now that we know it works, try to have fun.


## requirements.txt

To produce the list of needed packages use:

```
deps <- renv::dependencies()

pkgs <- unique(deps$Package)

writeLines(pkgs, "requirements.txt")
``` 

Ta-duh