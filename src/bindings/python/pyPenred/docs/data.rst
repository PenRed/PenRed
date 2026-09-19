Data Module
===========

.. automodule:: pyPenred.data
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: results, results1D, results2D, results3D, results4D, results5D,
                     results6D, results7D, results8D, results9D, results10D

Results classes
---------------

The module exposes the following classes with different maximum dimensionality:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Class
     - Stored data
   * - ``results1D``
     - 1-dimensional buffer
   * - ``results2D``
     - 2-dimensional buffer
   * - ``results3D``
     - 3-dimensional buffer
   * - ``results4D``
     - 4-dimensional buffer
   * - ``results5D``
     - 5-dimensional buffer
   * - ``results6D``
     - 6-dimensional buffer
   * - ``results7D``
     - 7-dimensional buffer
   * - ``results8D``
     - 8-dimensional buffer
   * - ``results9D``
     - 9-dimensional buffer
   * - ``results10D``
     - 10-dimensional buffer
   * - ``results``
     - buffer with maximum dimensions available (30)

They all share the same interface; only the dimensionality of the
underlying buffer differs. See :class:`~pyPenred.data.results` below
for the full method reference.

.. autoclass:: pyPenred.data.results
   :members:
   :undoc-members:
   :show-inheritance:

Dimension-specific methods
--------------------------

``results1D``
^^^^^^^^^^^^^

.. autoclass:: pyPenred.data.results1D
   :members: extractValue1D
   :undoc-members:

``results2D``
^^^^^^^^^^^^^

.. autoclass:: pyPenred.data.results2D
   :members: extractValue2D
   :undoc-members:

``results3D``
^^^^^^^^^^^^^

.. autoclass:: pyPenred.data.results3D
   :members: extractValue3D
   :undoc-members:
