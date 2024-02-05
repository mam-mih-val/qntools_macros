# QnTools-macros
## Common information

This project is an interface for QnTools framework for corrections and correlations in flow computation. The interace is set by root-macros using ROOT::RDataFrame. The examples can be found in the directory macro.

## Compilation

To compile this project the installed QnTools framework is required. We recomend using the ROOT version no older than 6.26. The QnTools framework can be obtained via the [link](https://github.com/FlowNICA/QnTools/tree/master).

When the QnTools is installed, export the framework
```
  $ export QnTools_DIR=/path/to/QnTools/install/lib/cmake/QnTools/
```
then proceed with instalation of the current project:
```
  $ cd /path/to/qntools_macros
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
```
