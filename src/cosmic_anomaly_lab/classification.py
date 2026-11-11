"""Configurable overlapping population labels and hypothesis applicability."""

from pathlib import Path

import yaml


class HypothesisRegistry:
    def __init__(self, path):
        self.config = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        self.populations = self.config["populations"]
        self.hypotheses = self.config["hypotheses"]

    def assess_labels(self, labels):
        """Labels are overlapping assessments, never a softmax partition."""
        unknown = set(labels) - set(self.populations)
        if unknown:
            raise ValueError(f"Unknown population labels: {unknown}")
        if any(not 0 <= confidence <= 1 for confidence in labels.values()):
            raise ValueError("Label confidences must lie in [0,1]")
        return {
            "labels": labels,
            "mutually_exclusive": False,
            "interpretation": "Confidence assessments are not proof of physical identity",
        }

    def applicability(self, hypothesis, available_observables):
        definition = self.hypotheses[hypothesis]
        missing = sorted(set(definition["requires"]) - set(available_observables))
        return {
            "hypothesis": hypothesis,
            "missing_observables": missing,
            "implemented_forward_model": definition["implemented_forward_model"],
            "testable": definition["implemented_forward_model"] and not missing,
        }
