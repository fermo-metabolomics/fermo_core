"""Runs the adduct annotation module.

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

import itertools
import logging
from typing import Self

from pydantic import BaseModel

from fermo_core.config.class_default_settings import DefaultMasses as Mass
from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Adduct,
    Annotations,
    Feature,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager
from fermo_core.utils.utility_method_manager import UtilityMethodManager

logger = logging.getLogger("fermo_core")


class AdductAnnotator(BaseModel):
    """Pydantic-based class to annotate "General Feature" objects with adduct info

    Attributes:
        params: ParameterManager object, holds user-provided parameters
        stats: Stats object, holds stats on molecular features and samples
        features: Repository object, holds "General Feature" objects
        samples: Repository object, holds "Sample" objects

    Sources:
        Blumer et al 2021 doi.org/10.1021/acs.jcim.1c00579
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository

    def return_features(self: Self) -> Repository:
        """Returns modified attributes from AdductAnnotator to the calling function

        Returns:
            Modified Feature Repository objects.
        """
        return self.features

    def run_analysis(self: Self):
        """Organizes calling of data analysis steps."""
        if self.params.PeaktableParameters.polarity == "positive":
            logger.info(
                "'AnnotationManager/AdductAnnotator': positive ion mode detected. "
                "Attempt to annotate for positive ion mode adducts."
            )
            for s_name in self.stats.samples:
                self.annotate_adducts_pos(s_name)
        else:
            logger.info(
                "'AnnotationManager/AdductAnnotator': negative ion mode detected. "
                "Attempt to annotate for negative ion mode adducts."
            )
            for s_name in self.stats.samples:
                self.annotate_adducts_neg(s_name)

        self.dereplicate_adducts()

    @staticmethod
    def add_adduct_info(feature: Feature) -> Feature:
        """Instantiate Annotation object instances

        Arguments:
            feature: Feature object to be modified

        Returns:
            (modified) Feature object
        """
        if feature.Annotations is None:
            feature.Annotations = Annotations()
        if feature.Annotations.adducts is None:
            feature.Annotations.adducts = []
        return feature

    def dereplicate_adducts(self: Self):
        """Combine identical adducts detected in different samples"""
        for f_id in self.stats.active_features:
            feature = self.features.get(f_id)
            if (
                feature.Annotations is not None
                and feature.Annotations.adducts is not None
                and len(feature.Annotations.adducts) > 0
            ):
                nonred_adducts = {}
                for adduct in feature.Annotations.adducts:
                    if adduct.partner_id not in nonred_adducts:
                        nonred_adducts[adduct.partner_id] = Adduct(
                            adduct_type=adduct.adduct_type,
                            partner_adduct=adduct.partner_adduct,
                            partner_id=adduct.partner_id,
                            partner_mz=adduct.partner_mz,
                            diff_ppm=adduct.diff_ppm,
                            sample_set={adduct.sample},
                        )
                    else:
                        nonred_adducts[adduct.partner_id].sample_set.add(adduct.sample)

                feature.Annotations.adducts = [val for _, val in nonred_adducts.items()]
                self.features.modify(f_id, feature)

    def annotate_adducts_neg(self: Self, s_name: str | int):
        """Pairwise compare features per sample, assign adducts info for negative mode

        Calculates overlap of features (peaks) by simplifying them to
        one-dimensional vectors. Consider two peaks A and B with A(start, stop)
        and B(start, stop). If any True in A_stop < B_start OR B_stop < A_start,
        peaks do NOT overlap. Base assumption is that one of the features is the
        [M-H]- ion.

        Arguments:
            s_name: a sample identifier

        """
        sample = self.samples.get(s_name)
        feature_set = sample.feature_ids.intersection(self.stats.active_features)

        if len(feature_set) == 0:
            logger.warning(
                f"'AnnotationManager/AdductAnnotator': no features to compare for "
                f"sample '{s_name}' - SKIP"
            )
            return

        f_pairs = itertools.combinations(feature_set, 2)
        for pair in f_pairs:
            feat1 = sample.features[pair[0]]
            feat2 = sample.features[pair[1]]

            if not (feat1.rt_stop < feat2.rt_start or feat2.rt_stop < feat1.rt_start):
                if self.chloride_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.chloride_adduct(feat1.f_id, feat2.f_id, s_name):
                    continue
                elif self.double_dimer_pair_neg(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_dimer_pair_neg(feat1.f_id, feat2.f_id, s_name):
                    continue
                elif self.bicarbonate_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.bicarbonate_adduct(feat1.f_id, feat2.f_id, s_name):
                    continue
                elif self.tfa_adduct(feat1.f_id, feat2.f_id, s_name) or self.tfa_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ):
                    continue
                elif self.acetate_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.acetate_adduct(feat1.f_id, feat2.f_id, s_name):
                    continue

    def annotate_adducts_pos(self: Self, s_name: str | int):
        """Pairwise compare features per sample, assign adducts info for positive mode

        Calculates overlap of features (peaks) by simplifying them to
        one-dimensional vectors. Consider two peaks A and B with A(start, stop)
        and B(start, stop). If any True in A_stop < B_start OR B_stop < A_start,
        peaks do NOT overlap. Base assumption is that one of the two features is
        the [M+H]+ adduct.

        Arguments:
            s_name: a sample identifier
        """
        sample = self.samples.get(s_name)
        feature_set = sample.feature_ids.intersection(self.stats.active_features)

        if len(feature_set) == 0:
            logger.warning(
                f"'AnnotationManager/AdductAnnotator': no features to compare for "
                f"sample '{s_name}' - SKIP"
            )
            return

        f_pairs = itertools.combinations(feature_set, 2)
        for pair in f_pairs:
            feat1 = sample.features[pair[0]]
            feat2 = sample.features[pair[1]]
            if not (feat1.rt_stop < feat2.rt_start or feat2.rt_stop < feat1.rt_start):
                if self.sodium_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.sodium_adduct(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.dimer_sodium_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.dimer_sodium_adduct(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.triple_h_adduct(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.triple_h_adduct(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.plus1_isotope(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.plus1_isotope(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.plus2_isotope(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.plus2_isotope(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.plus3_isotope(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.plus3_isotope(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.plus4_isotope(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.plus4_isotope(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.plus5_isotope(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.plus5_isotope(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.double_plus1(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_plus1(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.double_plus2(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_plus2(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.double_plus3(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_plus3(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.double_plus4(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_plus4(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.double_plus5(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.double_plus5(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.iron56(feat1.f_id, feat2.f_id, s_name) or self.iron56(
                    feat2.f_id, feat1.f_id, s_name
                ):
                    continue
                elif self.dimer_double(
                    feat1.f_id, feat2.f_id, s_name
                ) or self.dimer_double(feat2.f_id, feat1.f_id, s_name):
                    continue
                elif self.ammonium(feat1.f_id, feat2.f_id, s_name) or self.ammonium(
                    feat2.f_id, feat1.f_id, s_name
                ):
                    continue
                elif self.potassium(feat1.f_id, feat2.f_id, s_name) or self.potassium(
                    feat2.f_id, feat1.f_id, s_name
                ):
                    continue
                elif self.water_add(feat1.f_id, feat2.f_id, s_name) or self.water_add(
                    feat2.f_id, feat1.f_id, s_name
                ):
                    continue
                elif self.water_loss(feat1.f_id, feat2.f_id, s_name) or self.water_loss(
                    feat2.f_id, feat1.f_id, s_name
                ):
                    continue

    def sodium_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+Na]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz - Mass().H + Mass().Na), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+Na]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+Na]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def dimer_sodium_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [2M+Na]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((2 * (mh_ion.mz - Mass().H)) + Mass().Na), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[2M+Na]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[2M+Na]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def triple_h_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+3H]3+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + Mass().H) / 3), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+3H]3+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+3H]3+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def plus1_isotope(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+1+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + (1 * Mass().C13_12)), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+1+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+1+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def plus2_isotope(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+2+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + (2 * Mass().C13_12)), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+2+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+2+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def plus3_isotope(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+3+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + (3 * Mass().C13_12)), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+3+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+3+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def plus4_isotope(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+4+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + (4 * Mass().C13_12)), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+4+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+4+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def plus5_isotope(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+5+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + (5 * Mass().C13_12)), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+5+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+5+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_plus1(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+1+2H]2+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + (1 * Mass().C13_12)) / 2),
            adduct.mz,
            adduct.f_id,
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+1+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+1+2H]2+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_plus2(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+2+2H]2+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + (2 * Mass().C13_12)) / 2),
            adduct.mz,
            adduct.f_id,
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+2+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+2+2H]2+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_plus3(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+3+2H]2+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + (3 * Mass().C13_12)) / 2),
            adduct.mz,
            adduct.f_id,
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+3+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+3+2H]2+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_plus4(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+4+2H]2+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + (4 * Mass().C13_12)) / 2),
            adduct.mz,
            adduct.f_id,
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+4+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+4+2H]2+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_plus5(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+5+2H]2+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H + (5 * Mass().C13_12)) / 2),
            adduct.mz,
            adduct.f_id,
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+5+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+5+2H]2+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def iron56(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+56Fe-2H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz - (3 * Mass().H) + Mass().Fe56), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+56Fe-2H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+56Fe-2H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def dimer_double(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+2H]2+ and [2M+H]+ adducts, add information

        Consider two overlapping peaks A and B
        peak A with m/z 1648.47;
        peak B with m/z 824.74.
        If A is assumed [M+H]+, B would be [M+2H]2+
        If B is assumed [M+H]+, A would be [2M+H]+
        Thus, assignment is performed for [M+2H]2+ and [2M+H]+ in parallel,
        since M cannot be determined without isotopic data.

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((mh_ion.mz + Mass().H) / 2), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[2M+H]+",
                    partner_adduct="[M+2H]2+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+2H]2+",
                    partner_adduct="[2M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def ammonium(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+NH4]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz - Mass().H + Mass().NH4), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+NH4]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+NH4]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def potassium(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+K]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz - Mass().H + Mass().K), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+K]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+K]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def water_add(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+H2O+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz + Mass().H2O), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M+H2O+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H2O+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def water_loss(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M-H2O+H]+ adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        mh_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (mh_ion.mz - Mass().H2O), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(mh_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+H]+",
                    partner_adduct="[M-H2O+H]+",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-H2O+H]+",
                    partner_adduct="[M+H]+",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def chloride_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+Cl]- adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        m_h_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (m_h_ion.mz + Mass().H + Mass().Cl35), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(m_h_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-H]-",
                    partner_adduct="[M+Cl]-",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+Cl]-",
                    partner_adduct="[M-H]-",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def double_dimer_pair_neg(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M-2H]2- and [2M-H]- adduct pair, add information

        Consider two overlapping peaks A and B:
        -peak A with m/z 1648.47;
        -peak B with m/z 823.73.
        If A is assumed [M-H]-, B would be [M-2H]2-
        If B is assumed [M-H]-, A would be [2M-H]-
        Thus, assignment is performed for [M-2H]2- and [2M-H]- in parallel,
        since M cannot be determined without isotopic data.

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        m_h_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            ((m_h_ion.mz - Mass().H) / 2), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(m_h_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[2M-H]-",
                    partner_adduct="[M-2H]2-",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-2H]2-",
                    partner_adduct="[2M-H]-",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def bicarbonate_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+HCO2]- adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        m_h_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (m_h_ion.mz + Mass().H + Mass().HCO2), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(m_h_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-H]-",
                    partner_adduct="[M+HCO2]-",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+HCO2]-",
                    partner_adduct="[M-H]-",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def tfa_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+TFA-H]- (trifluoroacetate) adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        m_h_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (m_h_ion.mz + Mass().H + Mass().TFA), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(m_h_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-H]-",
                    partner_adduct="[M+TFA-H]-",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+TFA-H]-",
                    partner_adduct="[M-H]-",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False

    def acetate_adduct(self: Self, feat1: int, feat2: int, s_name: str) -> bool:
        """Determination of [M+HAc-H]- (acetate) adduct, add information

        Arguments:
            feat1: feature 1 identifier
            feat2: feature 2 identifier
            s_name: the sample identifier

        Returns:
            A bool indicating the outcome
        """
        m_h_ion = self.features.get(feat1)
        adduct = self.features.get(feat2)

        ppm = UtilityMethodManager.mass_deviation(
            (m_h_ion.mz + Mass().H + Mass().Ac), adduct.mz, adduct.f_id
        )
        if ppm < self.params.AdductAnnotationParameters.mass_dev_ppm:
            mh_ion = self.add_adduct_info(m_h_ion)
            adduct = self.add_adduct_info(adduct)
            mh_ion.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M-H]-",
                    partner_adduct="[M+HAc-H]-",
                    partner_id=adduct.f_id,
                    partner_mz=adduct.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            adduct.Annotations.adducts.append(
                Adduct(
                    adduct_type="[M+HAc-H]-",
                    partner_adduct="[M-H]-",
                    partner_id=mh_ion.f_id,
                    partner_mz=mh_ion.mz,
                    diff_ppm=ppm,
                    sample=s_name,
                )
            )
            self.features.modify(feat1, mh_ion)
            self.features.modify(feat2, adduct)
            return True
        else:
            return False
