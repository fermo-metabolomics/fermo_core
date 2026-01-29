"""Run the phenotype qualitative data assignment.

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
from statistics import mean, median
from typing import Optional, Self

import numpy as np
from pydantic import BaseModel
from scipy.stats import ranksums, ttest_ind

from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Annotations,
    Feature,
    Phenotype,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager

logger = logging.getLogger("fermo_core")


class PhenQualAssigner(BaseModel):
    """Pydantic-based class to run qualitative phenotype assignment

    Attributes:
        params: ParameterManager object, holds user-provided parameters
        stats: Stats object, holds stats on molecular features and samples
        features: Repository object, holds "General Feature" objects
        samples: Repository object holding Sample objects
        f_ids_intersect: the intersection of positive and negative feature IDs
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository
    f_ids_intersect: Optional[set] = None

    def return_values(self: Self) -> tuple[Stats, Repository]:
        """Return the modified objects to the calling function

        Returns:
            The modified Stats and Repository objects
        """
        return self.stats, self.features

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

    def collect_sets(self: Self):
        """Collect sets of active and inactive features and assign actives"""
        f_ids_all_actives = set()
        for s_id in set(self.stats.samples).difference(
            self.stats.phenotypes[0].s_negative
        ):
            sample = self.samples.get(s_id)
            f_ids_all_actives.update(sample.feature_ids)
        f_ids_all_actives.difference_update(self.stats.GroupMData.blank_f_ids)

        f_ids_all_inactives = set()
        for s_id in self.stats.phenotypes[0].s_negative:
            sample = self.samples.get(s_id)
            f_ids_all_inactives.update(sample.feature_ids)
        f_ids_all_inactives.difference_update(self.stats.GroupMData.blank_f_ids)

        f_ids_only_actives = f_ids_all_actives.difference(f_ids_all_inactives)
        self.stats.phenotypes[0].f_ids_positive.update(f_ids_only_actives)
        for f_id in f_ids_only_actives:
            feature = self.features.get(f_id)
            feature = self.add_annotation_attribute(feature=feature)
            feature.Annotations.phenotypes.append(
                Phenotype(
                    score=1.0, format="qualitative", descr="only in positive samples"
                )
            )
            self.features.modify(f_id, feature)

        self.f_ids_intersect = f_ids_all_actives.intersection(f_ids_all_inactives)

    def get_value(self: Self, f_id: int, sample_ids: set) -> list:
        """Retrieve values based on area or height

        Arguments:
            f_id: a feature ID
            sample_ids: the sample Ids to retrieve

        Returns:
            A list of collected values

        Raises:
            RuntimeError: empty list detected - no values determined
        """
        feature = self.features.get(f_id)
        values_s_ids = []

        if self.params.PhenoQualAssgnParams.value == "area":
            for entry in feature.area_per_sample:
                if entry.s_id in sample_ids:
                    values_s_ids.append(entry.value)
        else:
            for entry in feature.height_per_sample:
                if entry.s_id in sample_ids:
                    values_s_ids.append(entry.value)

        if len(values_s_ids) == 0:
            raise RuntimeError(
                "'PhenQualAssigner': empty list detected, would lead "
                "to ZeroDivisionError."
            )
        else:
            return values_s_ids

    def determine_active(
        self, feature: Feature, fct: float, act: list, inact: list
    ) -> Feature:
        """Determine active based on fold-change and test"""
        if fct < self.params.PhenoQualAssgnParams.factor:
            return feature

        phenotype = Phenotype(
            score=1.0,
            format="qualitative",
            descr=f"Fold-difference: '{round(fct, 2)}'",
        )

        if (
            len(act) > 3
            and len(inact) > 3
            and self.params.PhenoQualAssgnParams.p_val_cutoff != 0
        ):
            log_act = np.log10(act)
            log_inact = np.log10(inact)

            match self.params.PhenoQualAssgnParams.test:
                case "Welsh":
                    stat, pval = ttest_ind(log_act, log_inact, equal_var=False)
                    if pval < self.params.PhenoQualAssgnParams.p_val_cutoff:
                        phenotype.descr = f"Fold-difference: '{round(fct, 2)}', Welsh's t-test statistic: '{round(stat, 2)}'"
                        phenotype.p_value = pval
                case "Wilcoxon":
                    stat, pval = ranksums(log_act, log_inact)
                    if pval < self.params.PhenoQualAssgnParams.p_val_cutoff:
                        phenotype.descr = f"Fold-difference: '{round(fct, 2)}', Wilcoxon rank-sum test statistic: '{round(stat, 2)}'"
                        phenotype.p_value = pval
                case _:
                    logger.warning(
                        f"PhenQualAssigner: test '{self.params.PhenoQualAssgnParams.test}' not recognized - SKIP"
                    )

        feature = self.add_annotation_attribute(feature=feature)
        feature.Annotations.phenotypes.append(phenotype)
        self.stats.phenotypes[0].f_ids_positive.add(feature.f_id)
        return feature

    def bin_intersection(self: Self):
        """Bin the intersection between positive and negative f_ids based on factor

        Raises:
            RuntimeError: unsupported algorithm
        """
        for f_id in self.f_ids_intersect:
            feature = self.features.get(f_id)

            s_actv = set(self.stats.samples).difference(
                self.stats.phenotypes[0].s_negative
            )
            s_actv = s_actv.intersection(feature.samples)
            vals_act = self.get_value(f_id, s_actv)

            s_inactv = self.stats.phenotypes[0].s_negative
            s_inactv = s_inactv.intersection(feature.samples)
            vals_inact = self.get_value(f_id, s_inactv)

            match self.params.PhenoQualAssgnParams.algorithm:
                case "minmax":
                    feature = self.determine_active(
                        feature=feature,
                        fct=min(vals_act) / max(vals_inact),
                        act=vals_act,
                        inact=vals_inact,
                    )
                case "mean":
                    feature = self.determine_active(
                        feature=feature,
                        fct=mean(vals_act) / mean(vals_inact),
                        act=vals_act,
                        inact=vals_inact,
                    )
                case "median":
                    feature = self.determine_active(
                        feature=feature,
                        fct=median(vals_act) / median(vals_inact),
                        act=vals_act,
                        inact=vals_inact,
                    )
                case _:
                    raise RuntimeError("'PhenQualAssigner': Unsupported algorithm.")

            self.features.modify(f_id, feature)

    def run_analysis(self: Self):
        """Run the phenotype annotation analysis

        Raise:
            RuntimeError: self.stats.phenotypes not specified
        """
        if self.stats.phenotypes is None:
            raise RuntimeError(
                "'PhenQualAssigner': self.stats.phenotypes not specified - SKIP"
            )

        self.collect_sets()
        self.bin_intersection()
