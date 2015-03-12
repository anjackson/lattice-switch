# Introduction #

The code is not that user friendly, but the supplied [scripts](https://code.google.com/p/lattice-switch/source/browse/#svn%2Ftrunk%2Fscripts) and [example configuration files](https://code.google.com/p/lattice-switch/source/browse/#svn%2Ftrunk%2Fruns) should help.


# Basic invocation #

The code is invoked like this, with all supplied classes on the classpath:
```
java -Xmx1024m -cp $CLASSPATH net.anjackson.physics.ls.LSSim input.par
```

The [LSSim class](https://code.google.com/p/lattice-switch/source/browse/trunk/src/net/anjackson/physics/ls/LSSim.java) is the main controlled, and all simulation parameters are controlled through the supplied parameter file. Some examples may help:
  * [An example config file, for hard-sphere FCC-HCP simulation](http://lattice-switch.googlecode.com/svn/trunk/runs/hs-fcchcp/hsls.par).
  * [A template configuration file for NVT hard-sphere simulation](https://code.google.com/p/lattice-switch/source/browse/trunk/scripts/hsls.nvt.par.template).
  * [A template configuration file for NPT hard-sphere simulation](https://code.google.com/p/lattice-switch/source/browse/trunk/scripts/hsls.npt.par.template).
  * [A rather complicated example of a script that runs many simulations to explore the parameter space.](https://code.google.com/p/lattice-switch/source/browse/trunk/scripts/nvt-df-run.py)