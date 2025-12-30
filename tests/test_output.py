"""Tests for output formatting functionality."""

from __future__ import annotations

import io

from yallhap.classifier import HaplogroupCall, QCScores, SNPStats
from yallhap.cli import _write_result


class TestTSVOutput:
    """Tests for TSV output format."""

    def test_tsv_output_format(self) -> None:
        """TSV output has correct header and data columns."""
        result = HaplogroupCall(
            sample="TEST001",
            haplogroup="R-M343",
            confidence=0.95,
            qc_scores=QCScores(
                qc1_backbone=0.98,
                qc2_terminal=0.95,
                qc3_path=0.92,
                qc4_posterior=0.90,
            ),
            path=["ROOT", "R", "R-M343"],
            defining_snps=["M343"],
            alternatives={},
            snp_stats=SNPStats(
                informative_tested=100,
                derived=50,
                ancestral=40,
                missing=10,
            ),
            reference="grch38",
            tree_version="YFull",
        )

        output = io.StringIO()
        _write_result(result, output, "tsv")

        lines = output.getvalue().strip().split("\n")
        assert len(lines) == 2  # Header + data

        # Check header
        header = lines[0].split("\t")
        assert header == [
            "sample",
            "haplogroup",
            "confidence",
            "qc1",
            "qc2",
            "qc3",
            "qc4",
            "derived",
            "ancestral",
            "missing",
        ]

        # Check data
        data = lines[1].split("\t")
        assert data[0] == "TEST001"
        assert data[1] == "R-M343"
        assert data[2] == "0.9500"
        assert data[3] == "0.9800"  # qc1
        assert data[4] == "0.9500"  # qc2
        assert data[5] == "0.9200"  # qc3
        assert data[6] == "0.9000"  # qc4
        assert data[7] == "50"  # derived
        assert data[8] == "40"  # ancestral
        assert data[9] == "10"  # missing

    def test_json_output_format(self) -> None:
        """JSON output produces valid JSON."""
        import json

        result = HaplogroupCall(
            sample="TEST001",
            haplogroup="R-M343",
            confidence=0.95,
            qc_scores=QCScores(),
            path=["ROOT", "R-M343"],
            defining_snps=["M343"],
            alternatives={},
            snp_stats=SNPStats(derived=10, ancestral=5, missing=2),
            reference="grch38",
            tree_version="YFull",
        )

        output = io.StringIO()
        _write_result(result, output, "json")

        # Should be valid JSON
        data = json.loads(output.getvalue())
        assert data["sample"] == "TEST001"
        assert data["haplogroup"] == "R-M343"
        assert data["confidence"] == 0.95
