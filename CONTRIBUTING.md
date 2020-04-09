# Contributing to CosTuuM

Below is a set of guidelines for contributing to CosTuuM. By following these guidelines, you help us to tackle potential issues in an efficient way, or make sure pull requests can be accepted in a timely manner. Feel free to ignore these guidelines, or propose changes to this document by opening a pull request.

## Table of contents

 - [Code of conduct](#code-of-conduct)
 - [Direct support](#direct-support)
 - [How can I contribute?](#how-can-i-contribute)
 - [Guidelines for new issues](#guidelines-for-new-issues)
   - [Bug reports](#bug-reports)
   - [Requests for new features](#requests-for-new-features)
 - [Guidelines for code contributions](#guidelines-for-code-contributions)
   - [C++ code style](#c++-code-style)
   - [Python code style](#python-code-style)
   - [Style for other languages](#style-for-other-languages)
   - [General contribution guidelines](#general-contribution-guidelines)
 - [License and disclaimer](#license-and-disclaimer)

## Code of conduct

We aim to develop and maintain CosTuuM in an open and constructive environment, according to the rules set out in our [code of conduct](CODE_OF_CONDUCT.md). Please follow these rules and report unacceptable behaviour to us.

## Direct support

If you are unfamiliar with the way github works, or you want to use CosTuuM for your own project, but you don't have experience with C++/CMake/Pyton/any other technology on which CosTuuM depends, you are welcome to contact us directly on bert.vandenbroucke@ugent.be. We will try to help you as soon as possible.

## How can I contribute?

There are various ways to contribute to CosTuuM: by testing and providing us with useful bug reports, by proposing additional features, or by directly contributing new code.

To provide bug reports or to provide additional features, please open a new issue if possible. This allows us to keep track of things on github in a structured manner. Please try to adhere to the [issue guidelines](#guidelines-for-new-issues) below.

To contribute new code, please fork this repository, make the necessary changes, and then open a pull request. Make sure you follow the [rules for code contributions](#guidelines-for-code-contributions) below.

## Guidelines for new issues

There are two important types of issues which can be useful for us:
 - bug reports: issues that notify us of unwanted or unexpected code behaviour in the current stable version of the code
 - requests for new features: issues that relate to new features that are not part of the current code, including features that might improve the user friendliness of the code

It is important to clearly distinguish between the two, as bug reports evidently have a higher priority. To this end, you can use the labelling system: please mark bug reports as `bug`, and use any other label for new feature requests.

### Bug reports

In order to efficiently solve bugs, please try to provide us with as much information as possible to solve the bug:
 - a minimal example that causes the bug
 - information about the system you are running on
 
If you are unable to provide a minimal example (because e.g. your use case depends on a large input data file), please try to recreate the bug in a smaller example. If that does not work, please indicate why you are unable to provide a minimal example, and be prepared to provide us with additional data if requested.

If you encounter a bug, please reconfigure the code with `-DACTIVATE_ASSERTIONS=True` and rerun. This will activate a range of run time checks that might help track down the origin of the bug.

### Requests for new features

Try to be as detailed as possible about your request:
 - tell us why you need the new feature, so that we can assign a priority to its implementation
 - give us a clear example of how you would use the feature
 - optionally, provide us with your own ideas of how to implement the feature

We will consider your request as soon as possible, and get in touch with you. Please be prepared to provide us with additional information, e.g. example input.

## Guidelines for code contributions

Below are detailed guidelines for code contributions. New contributions should always be made to a fork of the original repository or in a seperate branch, and should be incorporated in the main code via a pull request to the `master` branch. All pull requests need to be approved by the repository administrator. Please make sure new code is up to date with the latest `master` branch by using `git rebase` (avoid using `git merge`) before opening the pull request, and run the code styling script (see below).

### C++ code style

In order to facilitate code comparisons across branches and commits, we make use of `clang-format`, an automated code styling tool. The default `clang-format` style file is part of the repository, and we provide a `bash` script that automatically applies this style to all C++ source code files. We kindly ask you to run this script before you commit any new code to the repository, and definitely when you open a new pull request. Note that this script uses `clang-format-6.0`, as there are small differences between different versions of the tool. The default style is not up for discussion, as most choices were made for good reasons.

In addition to the automatic styling, there are some additional rules:
 - all code blocks should be treated as code blocks, i.e. every `if`, `while` and `for` statement should be followed by a `{}` block, even if they consist of a single line of code.
 - single line conditional statements (`(condition) ? code for true : code for false`) are only allowed if they are truly single line statements or if they are necessary to initialize a `const` variable. In all other cases, a proper `if` statement is preferred.
 - `const` statements should be used wherever applicable. In fact, all variables should be declared as `const` unless they truly change. The same holds for class member routines that do not change the value of class variables.
 - global variables are not allowed in any part of the code, except in the code that is used to store compilation and configuration info.
 - all files, classes, member functions and member variables should be documented inline using `doxygen`. All of these should have at least a `@brief` statement giving a short description of that specific code unit. Documentation for member functions should always be given in the file that implements the function (so the `.cpp` file if both a `.hpp` and `.cpp` file are used for a class). All parameters and the return value of a member function should be fully documented; physical quantities (in parameters, return values and class members) should mention the units they use, with SI being favoured (unless not using SI units is more efficient).
 - `#ifdef` statements (and other preprocessor directives) should be avoided as much as possible. If their use is required, we advise using macros in the code body that are defined at the top of the file.
 - fast integer types (`int_fast32_t`...) should be used instead of default integer types, and should not have more precision than necessary. When looping over vectors or other standard library containers `size_t` should be used. `unsigned` integer types should be used for values that are strictly positive.
 - variable names should be clear and should use lowercase and underscores (no camel-case). Member variables should start with a an underscore (e.g. `_member_variable`). Template names should use upper case (e.g. `DATA_TYPE`). Class names should start with a capital and use camel-case; class names for deriving classes should extend the name of the parent class (e.g. `DraineDustProperties` as an implementation of a `DustProperties`). The same rules apply to `namespace` and `enum` names. Function names follow the normal variable styling. `enum` elements should use all caps, and should have names that start with the name of the corresponding `enum`.

### Python code style

Although the body of the CosTuuM library is and should be written in C++, Python can be used for unit testing and should be used for all benchmarks and example scripts. The same formatting script that automatically formats the C++ code will also automatically format Python scripts, using `black`. Again, we ask you to make sure all Python scripts have been formatted using this script before opening a pull request.

`black` formatting is a lot more restrictive than `clang-format` formatting, so generally it will impose a style that is consistent with what we consider to be good practice. Python naming conventions are the same as the corresponding C++ conventions. Some additional Python specific guidelines:
 - Python3 is the norm. Any old Python2 code should be converted using `2to3`.
 - modules or part of modules should be imported as a whole, e.g. `import numpy` and not `from numpy import *`.
 - imported modules can be aliased when convenient. We suggest the aliases `np` for `numpy`, `pl` for `matplotlib.pyplot` and `sp` for `scipy`.
 - scripts that read command line arguments should use `argparse` instead of providing a custom parsing routine. We generally encourage the development of general scripts that transform input into output, where both the input and output are command line arguments. This makes it easier to include these scripts in general workflows.
 - we strongly encourage the use of `numpy` and `scipy` to avoid explicit loops, e.g. for array operations or statistical analysis.

### Style for other languages

There is as of now no `clang-format` alternative for CMake or bash, so that we cannot enforce consistent styling for these languages. Please have a look at source code in those languages that is already part of the repository and try to use the same style. Most importantly:
 - use the same 80 character line limit in all files
 - try to mimick the `clang-format` styling for mathematical expressions: insert spaces in between variables and operands (`x + y` instead of `x+y` and same for `-`, `*` and `/`, but `x**y` and not `x ** y`) and split long lines in a logical and clear way, whenever possible.
 - use spaces instead of tabs and use 2 space indents where necessary
 - provide as much documentation as possible and at least as much as necessary

### General contribution guidelines

We aim to provide public, reproducible and stable code, which means adhering to the following rules:
 - all classes/functions should be unit tested. We use a CTest unit testing framework in which every class or namespace maps to a single unit test source file. Unit tests should cover as much code functionality as possible and should output clear error messages that help identify problems. Many examples are present in the `test` folder.
 - all mathematical expressions that are not immediately obvious should be properly documented in the `doxygen` documentation (`doxygen` has LaTeX support), and should cite relevant sources where applicable (algorithms based on Wikipedia articles should cite the relevant Wikipedia page; scientific publications are cited as e.g. Vandenbroucke et al. (2020) or Vandenbroucke & Wood (2018) - please make sure the reference points to a unique publication and insert the NASA ADS identifier or a DOI if possible).
 - proprietary code should not enter the public code if this is in conflict with the free software license of CosTuuM. If you are the owner of the proprietary code, you need to be aware that opening a pull request will make your code public. Please do not open a pull request if you do not want this.
 - new code should not duplicate existing functionality, i.e. if another class already implements the functionality you need, try to use that class or find a way to isolate that functionality in a new class that can be used by both your new class and the original class (it might be helpful if you open an issue to coordinate this)

## License and disclaimer

All code contributions will be public under the GNU Affero General Public License, and should therefore include the appropriate license banner (banners for various languages are included in a dedicated file in the repository). New code contributions should include appopriate author copyright information; code that results from suggested features will acknowledge the original author of the suggestion when possible.

You are free to use CosTuuM for your own work without contacting us, but we kindly ask you to cite the code paper (mentioned on the main repository page) whenever you do. We also suggest you include information about the specific version of the code that was used to any publication to improve reproducibility of the results.

Under the GNU Affero General Public License, no warranties can be given about the proper working of any of the code in this repository. So while we strive to ensure that CosTuuM produces accurate and reproducible results, we do not offer any guarantee that this is actually the case. We will try to offer support for any problems you might experience while running the code, but cannot guarantee that we will be able to do this in a timely and satisfactory manner. In the same way, we hope you can offer support for any code contributions you might make, but only if this suits you.
