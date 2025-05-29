"""Run the phenotype quantitative percentage data assignment.
Copyright (c) 2022 to present Mitja Maximilian Zdouc, PhD
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging
from typing import Self

import numpy as np
from pydantic import BaseModel
from scipy.stats import pearsonr, zscore
from statsmodels.stats.multitest import multipletests

from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Annotations,
    Feature,
    Phenotype,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats

logger = logging.getLogger("fermo_core")


class PhenQuantAss(BaseModel):
    """Pydantic-based class to run qualitative phenotype assignment

    Attributes:
        coeff_cutoff: the correlation coefficient cutoff
        p_val_cutoff: the corrected p-value cutoff
        fdr_corr: the type of FDR correction
        mode: the type of quantitative data
        stats: Stats object, holds stats on molecular features and samples
        features: Repository object, holds "General Feature" objects
        samples: Repository object holding Sample objects
        relevant_f_ids: features detected in > 3 samples
        assays: a dict of dicts containing correlations and p-values per feature for later assignment
    """

    coeff_cutoff: float
    p_val_cutoff: float
    fdr_corr: str
    mode: str
    stats: Stats
    features: Repository
    samples: Repository
    relevant_f_ids: set = set()
    assays: dict = {}

    def return_values(self: Self) -> tuple[Stats, Repository]:
        """Return the modified objects to the calling function

        Returns:
            The modified Stats and Repository objects
        """
        return self.stats, self.features

    def run_analysis(self: Self):
        """Run the phenotype annotation analysis

        Raise:
            RuntimeError: self.stats.phenotypes not assigned
        """
        if self.stats.phenotypes is None:
            raise RuntimeError(
                "'PhenQuantAssigner': self.stats.phenotypes not assigned - SKIP"
            )

        self.find_relevant_f_ids()
        self.calculate_correlation()
        self.correct_p_val()
        self.assign_results()

    @staticmethod
    def add_annotation_attribute(feature: Feature) -> Feature:
        """Add annotation attribute to feature if not existing

        Arguments:
            feature: the Feature object to modify

        Returns:
            The modified feature object
        """
        if feature.Annotations is None:
            feature.Annotations = Annotations()
        if feature.Annotations.phenotypes is None:
            feature.Annotations.phenotypes = []
        return feature

    def find_relevant_f_ids(self: Self):
        """Determines features detected in > 3 samples

        Raises:
            RuntimeError: No relevant features (present in >3 samples) detected.
        """
        for f_id in self.stats.active_features:
            feature = self.features.get(f_id)
            if len(feature.samples) > 3:
                self.relevant_f_ids.add(f_id)
            else:
                logger.debug(
                    f"'PhenQuantAssigner': feature id '{f_id}' only detected in "
                    f"'{len(feature.samples)}' samples: exclude from correlation "
                    f"analysis."
                )

        if len(self.relevant_f_ids) == 0:
            raise RuntimeError(
                "'PhenQuantAssigner': No relevant features (detected in > 3 "
                "samples) detected - SKIP."
            )
        else:
            logger.info(
                f"'PhenQuantAssigner':' Nr tested features: {len(self.relevant_f_ids)}"
            )

    @staticmethod
    def valid_constant(values: list) -> bool:
        """Check for uniform values

        Cannot perform Pearson corr on constant values

        Args:
            values: the measured values (areas, activities)

        Returns:
            A bool indicating outcome of validation
        """
        if np.isnan(values).any() or len(set(values)) == 1:
            return False
        else:
            return True

    @staticmethod
    def pearson_percentage(areas: list, activs: list) -> tuple[float, float]:
        """Calculate regular pearson coefficient

        Args:
            areas: the feature areas per sample
            activs: the measured activities per sample

        Returns:
            the pearson score and the p-value
        """
        return pearsonr(zscore(areas), zscore(activs))

    @staticmethod
    def pearson_concentration(areas: list, activs: list) -> tuple[float, float]:
        """Calculate regular pearson coefficient

        Args:
            areas: the feature areas per sample
            activs: the measured activities per sample

        Returns:
            the pearson score and the p-value
        """
        activs = [(1 / val) if val != 0 else 0 for val in activs]
        return pearsonr(zscore(areas), zscore(activs))

    def calculate_correlation(self: Self):
        """Collect data and prepare calculation

        Raises:
            KeyError: unsupported type of input data
        """
        for assay in self.stats.phenotypes:
            self.assays[assay.category] = {}

        for f_id in self.relevant_f_ids:
            feature = self.features.get(f_id)
            f_areas = {val.s_id: val.value for val in feature.area_per_sample}

            for assay in self.stats.phenotypes:
                areas = []
                activs = []
                for s_id in feature.samples:
                    for measure in assay.s_phen_data:
                        if s_id == measure.s_id:
                            areas.append(f_areas[s_id])
                            activs.append(measure.value)

                if len(areas) < 4 or len(activs) < 4:
                    continue

                if not self.valid_constant(areas) or not self.valid_constant(activs):
                    continue

                if self.mode == "percentage":
                    pearson_s, p_val = self.pearson_percentage(
                        areas=areas, activs=activs
                    )
                elif self.mode == "concentration":
                    pearson_s, p_val = self.pearson_concentration(
                        areas=areas, activs=activs
                    )
                else:
                    raise KeyError(
                        f"Unsupported type of input data '{self.mode}' - SKIP"
                    )

                self.assays[assay.category][f_id] = {}
                self.assays[assay.category][f_id]["corr"] = pearson_s
                self.assays[assay.category][f_id]["raw_p"] = p_val
                self.assays[assay.category][f_id]["datatype"] = assay.datatype

    def correct_p_val(self) -> None:
        """Correct the calculated p-values with a user-defined algorithm"""
        for assay in self.assays:
            raw_p_values = [v["raw_p"] for v in self.assays[assay].values()]
            adjusted_p_values = multipletests(raw_p_values, method=self.fdr_corr)[1]

            for feature, adj_p in zip(self.assays[assay].keys(), adjusted_p_values):
                self.assays[assay][feature]["adjusted_p"] = adj_p

    def assign_results(self) -> None:
        """Assign the results to the feature repository"""

        def _check_permitted(score: float, p_corr: float) -> bool:
            if self.coeff_cutoff == 0 or self.p_val_cutoff == 0:
                return True
            elif score > self.coeff_cutoff and p_corr < self.p_val_cutoff:
                return True
            else:
                return False

        for f_id in self.relevant_f_ids:
            feature = self.features.get(f_id)
            feature = self.add_annotation_attribute(feature=feature)

            for num, assay in enumerate(self.assays):
                try:
                    data = self.assays[assay][f_id]
                    if _check_permitted(data["corr"], data["adjusted_p"]):
                        feature.Annotations.phenotypes.append(
                            Phenotype(
                                format=data["datatype"],
                                category=assay,
                                score=data["corr"],
                                p_value=data["raw_p"],
                                p_value_corr=data["adjusted_p"],
                            )
                        )
                        self.stats.phenotypes[num].f_ids_positive.add(f_id)
                except Exception as e:
                    logger.warning(
                        f"'PhenQuantAssigner': assignment for {e!s} failed: insufficient values to calculate correlation"
                    )

            self.features.modify(f_id, feature)
