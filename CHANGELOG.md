# Unreleased

# v1.7.0
## Features
- There is now an option to use a policy-based algorithm for the forward phase (also known as control maps). See `ForwardMode` and the Theory guide for more information.

## Fixes
- Fixed a bug breaking the split model mode when the unfeasibility flag is left empty.

# v1.6.5
## Features
- Added a check on the number of state variables returned by the model function.

## Fixes
- Fixed a check on the number of terminal state constraints.

# v1.6.4
## Features
- There is now a 'Display' option that can be used to adjust verbosity (suppress warnings and/or the progress bar)

# v1.6.3
## Features
- A warning is now issued if the system and cost function returns full unfeasibility at a given stage.
- There is a new property called totalCost containing the total cost (shock). 
- The code now only issues warnings if the optimization fails instead of errors.
- If the optimization fails, totalCost is set to infinity.

# v1.6.2
## Docs
- Expand and clarified the docs for the cart example.

## Fixes
- Fixed a code snippet in the hev_example which erroneously showed code from another example.
- Fixed a code snippet in the "Set up a basic problem" section where the final state constraints are described.
- Fixed a sample code snippet in the additional inputs documentation when not using exogenous inputs.

# v1.6.1

## Fixes
- Fixed a bug that could prevent the 'EnforceStateGrid' option from effectively enforcing a constraint on the state grid bounds.

# v1.6

## Features
- There is now an option to specify a custom terminal cost. An example for its usage has also been added.
- There is a new example illustrating the use of the level set method. 

## Fixes
- Fixed a bug that occasionally occured when using exogenous inputs in conjunction with safe mode.

## Misc
- I refactored most of the code to further modularize DynaProg's methods. This all happened "under the hood" and it should be useful for developers, while being irrelevant to regular users.

# v1.5.2

## Fixes
- Fixed a bug where the forward phase would occasionally fail when using safe mode and exogenous inputs.

# v1.5.1

## Fixes
- Fixed a bug where the forward simulation progress counter would stop if warnings arose.


# v1.5

## Features
- There is now an option to also store the value functions constructed in the DP backward phase. Check out `StoreValueFunction` in the Syntax guide for more info. 
- There are some new options to change the way the terminal state constraints are enforced. Check out `VFInitialization`, `VFFfactors` and `myInf` in the Syntax guide for more info.
- There is now an option to automatically enforce a constraint on the state variables so that they do not exceed the state grids. Note that this is now the default behavior. Check out and `EnforceStateGrid` in the Syntax guide for more info. 
- Online docs are now available at [https://fmiretti.github.io/DynaProg/](https://fmiretti.github.io/DynaProg/). 
- DynaProg objects now have a tag with the version number they were created with.
- Added more checks and warnings.

## Fixes
- Fixed a mismatch between the cart example and its documentation.
- Fixed an error occurring in the cumulative cost plot.
- Other minor fixes.