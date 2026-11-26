#!/usr/bin/env python3
"""Smoke tests for regular- and variable-grid pyPenred data classes.

Run with the pyPenred build/install directory on PYTHONPATH:

    python3 test_data_grids.py
"""

import os
import tempfile
import unittest

import numpy as np

from pyPenred import data


class ResultsGridTests(unittest.TestCase):
    def test_registered_dimensions(self):
        for dimension in range(1, 11):
            regular = getattr(data, "results{}D".format(dimension))()
            variable = getattr(data, "variableResults{}D".format(dimension))()
            self.assertEqual(regular.maxDim(), dimension)
            self.assertEqual(variable.maxDim(), dimension)

        self.assertEqual(data.results().maxDim(), 30)
        self.assertEqual(data.variableResults().maxDim(), 30)

    def test_regular_1d_fill_data_interpolation_and_io(self):
        values = np.array([1.0, 2.0, 4.0])
        sigma = np.full(values.shape, 0.1)
        result = data.results1D()
        result.fill(values, sigma, [(0.0, 3.0, "x")], "value")

        self.assertEqual(result.bins(0), 3)
        self.assertEqual(result.bins(), 3)
        self.assertEqual(result.title(), "")
        self.assertEqual(result.maxDim(), 1)

        read_values, read_sigma, info, value_header, description = result.data()
        np.testing.assert_array_equal(read_values, values)
        np.testing.assert_array_equal(read_sigma, sigma)
        self.assertEqual(tuple(info), (0.0, 3.0, "x"))
        self.assertEqual(value_header, "value")
        self.assertEqual(description, "")

        self.assertEqual(len(result.extractValue([1.5])), 2)
        self.assertEqual(len(result.extractValue1D(1.5)), 2)

        with tempfile.TemporaryDirectory() as directory:
            filename = os.path.join(directory, "regular.dat")
            result.save(filename)
            loaded = data.results1D()
            loaded.load(filename)
            loaded_values, loaded_sigma, loaded_info, _, _ = loaded.data()
            np.testing.assert_array_equal(loaded_values, values)
            np.testing.assert_array_equal(loaded_sigma, sigma)
            self.assertEqual(tuple(loaded_info), (0.0, 3.0, "x"))

    def test_variable_1d_edges_and_interpolation(self):
        values = np.array([1.0, 3.0, 5.0])
        sigma = np.full(values.shape, 0.2)
        edges = [0.0, 0.5, 2.0, 5.0]
        result = data.variableResults1D()
        result.fill(values, sigma, [(edges, "energy")], "counts")

        self.assertEqual(result.bins(0), 3)
        read_values, read_sigma, info, value_header, _ = result.data()
        np.testing.assert_array_equal(read_values, values)
        np.testing.assert_array_equal(read_sigma, sigma)
        np.testing.assert_array_equal(info[0], edges)
        self.assertEqual(info[1], "energy")
        self.assertEqual(value_header, "counts")
        self.assertEqual(len(result.extractValue([1.0])), 2)
        self.assertEqual(len(result.extractValue1D(1.0)), 2)

    def test_2d_operations_spectrum_and_variable_grid(self):
        values = np.arange(6.0).reshape(2, 3)
        sigma = np.full(values.shape, 0.1)
        regular = data.results2D()
        regular.fill(
            values,
            sigma,
            [(0.0, 2.0, "y"), (0.0, 3.0, "x")],
            "value",
        )

        self.assertEqual(regular.bins(0), 3)
        self.assertEqual(regular.bins(1), 2)
        self.assertEqual(len(regular.extractValue([1.0, 1.0])), 2)
        self.assertEqual(len(regular.extractValue2D(1.0, 1.0)), 2)

        spectrum, spectrum_sigma = regular.extractSpectrum1D(0, [1.0])
        self.assertEqual(spectrum.shape, (3,))
        self.assertEqual(spectrum_sigma.shape, (3,))

        difference = regular.subtract(regular)
        np.testing.assert_array_equal(difference.data()[0], np.zeros_like(values))
        ratio = regular.relative(regular)
        expected_ratio = np.divide(
            values,
            values,
            out=np.zeros_like(values),
            where=values != 0,
        )
        np.testing.assert_allclose(ratio.data()[0], expected_ratio)

        normalized = data.results1D()
        normalized.fill(
            np.array([1.0, 2.0]),
            np.array([0.1, 0.1]),
            [(0.0, 2.0, "x")],
        )
        normalized.normalize()
        np.testing.assert_allclose(normalized.data()[0].sum(), 1.0)

        variable = data.variableResults2D()
        variable.fill(
            values,
            sigma,
            [
                ([0.0, 1.0, 4.0], "y"),
                ([0.0, 0.2, 1.0, 3.0], "x"),
            ],
        )
        variable_values, _, variable_y, variable_x, _, _ = variable.data()
        np.testing.assert_array_equal(variable_values, values)
        np.testing.assert_array_equal(variable_y[0], [0.0, 1.0, 4.0])
        self.assertEqual(variable_y[1], "y")
        np.testing.assert_array_equal(variable_x[0], [0.0, 0.2, 1.0, 3.0])
        self.assertEqual(variable_x[1], "x")
        self.assertEqual(len(variable.extractValue2D(0.5, 0.5)), 2)

    def test_3d_dimension_conversions(self):
        values_2d = np.ones((2, 2, 1))
        result_2d = data.variableResults3D()
        result_2d.fill(
            values_2d,
            np.ones_like(values_2d),
            [
                ([0.0, 1.0, 2.0], "z"),
                ([0.0, 1.0, 2.0], "y"),
                ([0.0, 1.0], "x"),
            ],
        )
        self.assertIsInstance(result_2d.to2D(), data.variableResults2D)

        values_1d = np.ones((2, 1, 1))
        result_1d = data.results3D()
        result_1d.fill(
            values_1d,
            np.ones_like(values_1d),
            [
                (0.0, 1.0, "z"),
                (0.0, 1.0, "y"),
                (0.0, 2.0, "x"),
            ],
        )
        self.assertIsInstance(result_1d.to1D(), data.results1D)


class SplineTests(unittest.TestCase):
    def test_results_create_cubic_spline(self):
        result = data.results1D()
        result.fill(
            np.array([0.0, 1.0, 8.0, 27.0]),
            np.ones(4),
            [(0.0, 3.0, "x")],
        )

        spline = result.interpolate1D(0)
        self.assertIsInstance(spline, data.CubicSpline)
        self.assertAlmostEqual(spline.evaluate(spline.readX()[2]), 8.0)

        limited = result.interpolate1D(0, [(0, 0, 4)])
        self.assertIsInstance(limited, data.CubicSpline)
        self.assertAlmostEqual(limited.evaluate(limited.readX()[2]), 8.0)

    def test_results_create_bicubic_spline(self):
        values = np.fromfunction(
            lambda row, column: row + column,
            (4, 4),
            dtype=float,
        )
        result = data.results2D()
        result.fill(
            values,
            np.ones_like(values),
            [(0.0, 3.0, "y"), (0.0, 3.0, "x")],
        )

        spline = result.interpolate2D(0, 1)
        self.assertIsInstance(spline, data.BicubicSpline)
        self.assertAlmostEqual(spline.evaluate(0.0, 0.0), 0.0)

    def test_cubic_spline_scalar_and_bulk_evaluation(self):
        spline = data.CubicSpline()
        spline.init(
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 1.0, 8.0, 27.0],
            s1=0.0,
            sn=18.0,
        )

        self.assertTrue(spline.initialized())
        self.assertEqual(list(spline.readX()), [0.0, 1.0, 2.0, 3.0])
        self.assertAlmostEqual(spline.evaluate(1.5), 3.375)
        np.testing.assert_allclose(
            spline.evaluate(np.array([0.5, 1.5, 2.5])),
            [0.125, 3.375, 15.625],
        )

    def test_cubic_spline_rejects_invalid_input(self):
        spline = data.CubicSpline()
        with self.assertRaises(ValueError):
            spline.init([0.0, 1.0, 1.0, 3.0], [0.0, 1.0, 2.0, 3.0])

    def test_bicubic_spline_scalar_and_grid_evaluation(self):
        spline = data.BicubicSpline()
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 1.0, 2.0, 3.0]
        values = [x_value + y_value for y_value in y for x_value in x]
        spline.init(x, y, values)

        self.assertAlmostEqual(spline.evaluate(1.5, 2.5), 4.0)
        np.testing.assert_allclose(
            spline.evaluateGrid([0.5, 1.5], [0.5, 1.5]),
            [[1.0, 2.0], [2.0, 3.0]],
        )

    def test_tensor_product_spline_axis_order_and_grid_shape(self):
        axes = [
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 1.0, 2.0, 3.0],
            [0.0, 1.0, 2.0, 3.0],
        ]
        values = np.empty((4, 4, 4))
        for i2, z in enumerate(axes[2]):
            for i1, y in enumerate(axes[1]):
                for i0, x in enumerate(axes[0]):
                    values[i2, i1, i0] = x + 2.0 * y + 3.0 * z

        spline = data.TensorProductSpline()
        spline.init(axes, values.ravel())

        self.assertAlmostEqual(spline.evaluate([0.5, 1.5, 2.5]), 11.0)
        grid = spline.evaluateGrid(
            [[0.5, 1.5], [0.5, 1.5], [0.5, 1.5]]
        )
        self.assertEqual(grid.shape, (2, 2, 2))
        self.assertAlmostEqual(grid[0, 0, 0], 3.0)
        self.assertAlmostEqual(grid[0, 0, 1], 4.0)

    def test_results_create_tensor_product_spline(self):
        values = np.empty((4, 4, 4))
        for i2 in range(4):
            for i1 in range(4):
                for i0 in range(4):
                    values[i2, i1, i0] = i0 + 2.0 * i1 + 3.0 * i2

        result = data.results3D()
        result.fill(
            values,
            np.ones_like(values),
            [(0.0, 4.0, "z"), (0.0, 4.0, "y"), (0.0, 4.0, "x")],
        )

        spline = result.interpolateND()
        self.assertIsInstance(spline, data.TensorProductSpline)
        self.assertAlmostEqual(spline.evaluate([0.5, 1.5, 2.5]), 11.0)


if __name__ == "__main__":
    unittest.main()
