# wadeTools
Package version of my toolbox

## Installing the package
```
devtools::install_github("rogerswt/wadeTools", build_vignettes = TRUE)
```
## Purpose of this package
The ***wadeTools*** package is a collection of a wide variety of R functions that I've written over the years.  They accomplish big or small tasks that I find that I do over and over again, so as in any software development work, it makes sense to write it once (carefully and with an eye towards re-useability) and use it often.  

I had originally created a 'tools' directory which had task-specific subdirectories (e.g. visualization, compensation, transformation, etc.).  Any time I found I had a function that seemed to have broader utility than just the project at hand, I dumped it into the tools hierarchy.  Over time functions and categories accumulated.  Mostly these had to do specifically with flow cytometry data processing, but a few didn't.

I've only recently turned this collection into a formal R package.  Why?  Because I wanted to share these tools with **Cytomics Workshop** participants - people who are interested in acquiring or honing  skills in developing and applying advanced computational analysis methods for flow cytometry.  In earlier iterations of the workshop I simply shared the tools directory.  That's akin to dumping a box of random tools on the workbench and asking you to build a barn, with no *a priori* knowledge of what a plane, or a saw, or a chisel was used for, or how to use them.  A package has several monumental advantages over this *ad hoc* approach to sharing code:

1.  it's easy to distribute via Github,
1.  it installs and loads just like any other R package,
1.  each function is documented with a help page that shows its arguments and usage, and in many cases examples of its use,
1.  it comes with this vignette - essentially a cook's tour of the package and how to get started with it.