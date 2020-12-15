=========
Changelog
=========

Version 0.1.5
=============

what's new:
-----------

- Pre-release, heading towards V0.2.x

- Restructuring for package release

- Lots more documentation (still incomplete)

- Fixes to package configuration

Known problems
--------------

- Incomplete documentation, especially outside the ``yapcad.geom`` module.

- Occasional problems with complex boolean operations

- Incomplete functionality around 3D modeling

- Inconsistent inclusion of licensing boilerplate
  
Version 0.2.0
=============

what's new:
-----------

- First announced version of **yapCAD**. Yay!

- Added new ``boxcut`` example, showing a fully worked (if simple)
  parametric design system.

- Additional documentation updates and minor bugfixes.

Known problems
--------------

- Our `yapCAD readthedocs`_ documentation is missing the expanded
  documentation from submodules, which is a problem since much of
  **yapCAD**'s documentation is in the form of docstrings in the
  source.  I'm working on getting this sorted out.  In the mean time,
  you may want to build a local copy of the documentation as described
  in the main ``README`` file.   Or, checkout and read the source. 
  
- Incomplete documentation, especially outside the ``yapcad.geom`` module.

- Occasional problems with complex boolean operations.  A bug in the
  ``intersectXY`` method of the ``Boolean`` class.

- Incomplete functionality around 3D modeling

- Inconsistent inclusion of licensing boilerplate, other minor
  formatting issues.

.. _yapCAD readthedocs: https://yapcad.readthedocs.io/en/latest/index.html
