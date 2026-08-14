from __future__ import annotations

import builtins
import math
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from lcms_curation.qc import _spearman_rho, compute_abundance_metrics


class QCStatisticsValidation(unittest.TestCase):
    def test_spearman_monotonic_and_inverse(self) -> None:
        x = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(_spearman_rho(x, np.array([10.0, 20.0, 30.0, 40.0])), 1.0, places=12)
        self.assertAlmostEqual(_spearman_rho(x, np.array([40.0, 30.0, 20.0, 10.0])), -1.0, places=12)

    def test_spearman_uses_average_tie_ranks_and_pairwise_finite_values(self) -> None:
        self.assertAlmostEqual(
            _spearman_rho(np.array([1.0, 2.0, 2.0, 4.0]), np.array([1.0, 2.0, 3.0, 4.0])),
            0.9486832980505138,
            places=12,
        )
        self.assertAlmostEqual(
            _spearman_rho(np.array([1.0, 2.0, np.nan, 4.0]), np.array([10.0, 20.0, 30.0, 40.0])),
            1.0,
            places=12,
        )
        self.assertTrue(math.isnan(_spearman_rho(np.ones(4), np.arange(4.0))))

    def test_qc_drift_path_never_imports_scipy(self) -> None:
        frame = pd.DataFrame(
            {
                "Feature ID": ["feature-1"],
                "q1": [100.0],
                "q2": [110.0],
                "q3": [125.0],
                "q4": [140.0],
            }
        )
        metadata = pd.DataFrame(
            {
                "Include": [True] * 4,
                "Sample role": ["Pooled QC"] * 4,
                "Normalized abundance column": ["q1", "q2", "q3", "q4"],
                "Injection order": [1, 2, 3, 4],
                "Dilution factor": [1.0] * 4,
                "Technical replicate": [False] * 4,
                "Subject / biological unit": ["QC"] * 4,
            }
        )
        original_import = builtins.__import__

        def reject_scipy(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy" or name.startswith("scipy."):
                raise ModuleNotFoundError("SciPy deliberately unavailable during regression test")
            return original_import(name, globals, locals, fromlist, level)

        with patch("builtins.__import__", side_effect=reject_scipy):
            result = compute_abundance_metrics(frame, metadata)
        self.assertAlmostEqual(float(result.loc[0, "QC injection-order Spearman rho"]), 1.0, places=12)


if __name__ == "__main__":
    unittest.main()
