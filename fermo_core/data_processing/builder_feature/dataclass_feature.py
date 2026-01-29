"""Organize molecular feature data. Product of Feature Builder class.

Copyright (c) 2022-2023 Mitja Maximilian Zdouc, PhD

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
from typing import Any, Optional, Self

from pydantic import BaseModel

logger = logging.getLogger("fermo_core")


class Adduct(BaseModel):
    """A Pydantic-based class to represent adduct annotation information

    Attributes:
        adduct_type: describes type of adduct
        partner_adduct: adduct of partner feature
        partner_id: partner feature id based on which the adduct was found
        partner_mz: partner feature mz
        diff_ppm: the difference in ppm between the two features
        sample: the sample identifier
        sample_set: a set of samples in which Adduct was observed
    """

    adduct_type: str
    partner_adduct: str
    partner_id: int
    partner_mz: float
    diff_ppm: float
    sample: Optional[str] = None
    sample_set: Optional[set] = None

    def to_json(self: Self) -> dict:
        return {
            "adduct_type": self.adduct_type,
            "partner_adduct": self.partner_adduct,
            "partner_id": self.partner_id,
            "partner_mz": self.partner_mz,
            "diff_ppm": round(self.diff_ppm, 2),
            "samples": list(self.sample_set) if self.sample_set is not None else [],
        }


class Match(BaseModel):
    """A Pydantic-based class to represent library matching information

    Attributes:
        id: identifier of matched molecule
        library: name of library file matched molecule is from
        algorithm: algorithm used for matching
        score: score that led to a matching between feature and matched molecule
        mz: m/z ratio of matched molecule
        diff_mz: difference between m/z ratios of feature and matched molecule
        module: responsible for matching operation
        smiles: optional smiles string (ms2query)
        inchikey: optional inchi key (ms2query)
        npc_class: NPClassifier class of analog (ms2query)
    """

    id: Any
    library: str
    algorithm: str
    score: float
    mz: float
    diff_mz: float
    module: str
    smiles: Optional[str] = None
    inchikey: Optional[str] = None
    npc_class: Optional[str] = None

    def to_json(self: Self) -> dict:
        return {
            "id": self.id,
            "library": self.library,
            "algorithm": self.algorithm,
            "score": self.score,
            "mz": self.mz,
            "diff_mz": self.diff_mz,
            "module": self.module,
            "smiles": self.smiles if self.smiles is not None else "N/A",
            "inchikey": self.inchikey if self.inchikey is not None else "N/A",
            "npc_class": self.npc_class if self.npc_class is not None else "N/A",
        }


class NeutralLoss(BaseModel):
    """A Pydantic-based class to represent detected neutral losses in MS2 spectrum

    Attributes:
        id: neutral loss identifier
        loss_det: detected m/z
        loss_ex: expected m/z
        mz_frag: the fragment m/z corresponding to neutral loss (parent m/z - loss)
        diff: difference in ppm
    """

    id: str
    loss_det: float
    loss_ex: float
    mz_frag: float
    diff: float

    def to_json(self: Self) -> dict:
        return {
            "id": self.id,
            "det_loss": round(self.loss_det, 4),
            "exp_loss": self.loss_ex,
            "mz_frag": round(self.mz_frag, 4),
            "diff_ppm": round(self.diff, 2),
        }


class CharFrag(BaseModel):
    """A Pydantic-based class to represent characteristic ion fragments in MS2 spectrum

    Attributes:
        id: fragment identifier
        frag_det: detected m/z
        frag_ex: expected m/z
        diff: difference in ppm
    """

    id: str
    frag_det: float
    frag_ex: float
    diff: float

    def to_json(self: Self) -> dict:
        return {
            "id": self.id,
            "frag_det": round(self.frag_det, 4),
            "frag_ex": round(self.frag_ex, 4),
            "diff_ppm": round(self.diff, 2),
        }


class Phenotype(BaseModel):
    """A Pydantic-based class to represent phenotype information

    Attributes:
        format: the format of the phenotype file
        category: the assay category (column) if applicable
        descr: additional data if applicable
        score: the score calculated
        p_value: the calculated p-value if applicable
        p_value_corr: the corrected p-value if applicable

    """

    format: str
    category: Optional[str] = None
    descr: Optional[str] = None
    score: float
    p_value: Optional[float] = None
    p_value_corr: Optional[float] = None

    def to_json(self: Self) -> dict:
        return {
            "format": self.format,
            "category": self.category if self.category is not None else "N/A",
            "descr": self.descr if self.descr is not None else "N/A",
            "score": round(self.score, 6),
            "p_value": round(self.p_value, 10) if self.p_value is not None else "N/A",
            "p_value_corr": (
                round(self.p_value_corr, 10) if self.p_value_corr is not None else "N/A"
            ),
        }


class Annotations(BaseModel):
    """A Pydantic-based class to represent annotation information

    Attributes:
        adducts: list of Adduct objects representing putative adducts of this feature
        matches: list of Match objects repr. putative library matching hits
        losses: list of NeutralLoss objects annotating functional groups of feature
        fragments: list of CharFrag objects annotating characteristic ion fragments
        phenotypes: list of Phenotype objects if feature phenotype-associated
    """

    adducts: Optional[list] = None
    matches: Optional[list] = None
    losses: Optional[list] = None
    fragments: Optional[list] = None
    phenotypes: Optional[list] = None

    def sort_entries(self: Self, attr: str, score: str, direction: bool):
        """Sort the entries in 'attr' based on 'score' in descending order

        Arguments:
            attr: the attribute to target
            score: the score to sort for
            direction: False to sort low to high, True for reverse

        """
        if getattr(self, attr) is not None:
            setattr(
                self,
                attr,
                sorted(
                    getattr(self, attr),
                    key=lambda x: getattr(x, score),
                    reverse=direction,
                ),
            )

    def to_json(self: Self) -> dict:
        self.sort_entries(attr="adducts", score="diff_ppm", direction=False)
        self.sort_entries(attr="matches", score="score", direction=True)
        self.sort_entries(attr="losses", score="diff", direction=False)
        self.sort_entries(attr="fragments", score="diff", direction=False)
        self.sort_entries(attr="phenotypes", score="score", direction=True)
        return {
            "adducts": (
                [adduct.to_json() for adduct in self.adducts]
                if self.adducts is not None
                else []
            ),
            "matches": (
                [adduct.to_json() for adduct in self.matches]
                if self.matches is not None
                else []
            ),
            "losses": (
                [adduct.to_json() for adduct in self.losses]
                if self.losses is not None
                else []
            ),
            "fragments": (
                [adduct.to_json() for adduct in self.fragments]
                if self.fragments is not None
                else []
            ),
            "phenotypes": (
                [phenotype.to_json() for phenotype in self.phenotypes]
                if self.phenotypes is not None
                else []
            ),
        }


class SimNetworks(BaseModel):
    """A Pydantic-based class to represent spectral similarity network information

    Attributes:
        algorithm: name of algorithm
        network_id: an integer indicating the network feature is associated with
    """

    algorithm: str
    network_id: int

    def to_json(self: Self) -> dict:
        return {
            "algorithm": self.algorithm,
            "network_id": self.network_id,
        }


class SampleInfo(BaseModel):
    """A Pydantic-based class to represent sample->value information

    Attributes:
        s_id: identifier of sample
        value: indicating the respective value (area or height)
    """

    s_id: str
    value: float

    def to_json(self: Self) -> dict:
        return {
            "s_id": str(self.s_id),
            "value": float(self.value),
        }


class GroupFactor(BaseModel):
    """A Pydantic-based class to represent group factor (fold difference) information

    Attributes:
        group1: comparison group 1
        group2: comparison group 2
        factor: group factor (fold difference), with the greater number taken
    """

    group1: str
    group2: str
    factor: float

    def to_json(self: Self) -> dict:
        return {
            "group1": str(self.group1),
            "group2": str(self.group2),
            "factor": round(self.factor, 2),
        }


class Scores(BaseModel):
    """A Pydantic-based class to represent feature score information

    Attributes:
        phenotype: the highest phenotype correlation score (if any) across all assays
        novelty: putative novelty of the feature (compared against external data)
    """

    phenotype: Optional[float] = None
    novelty: Optional[float] = None

    def to_json(self) -> dict:
        return {
            "phenotype": (
                round(self.phenotype, 2) if self.phenotype is not None else 0.0
            ),
            "novelty": round(self.novelty, 2) if self.novelty is not None else 1.0,
        }


class Feature(BaseModel):
    """A Pydantic-based class to represent a molecular feature.

    Attributes:
        f_id: the integer ID of the molecular feature.
        mz: the precursor mass to charge ratio (m/z).
        rt: the retention time at peak apex.
        rt_start: the start of the peak in minutes.
        rt_stop: the stop of the peak in minutes.
        rt_range: the length of a molecular feature peak in minutes.
        trace_rt: the rt data points of the pseudo-chromatogram trace
        trace_int: the relative intensity data points of the pseudo-chromatogram trace
        fwhm: the feature width at half maximum intensity (peak width).
        intensity: the maximum intensity.
        rel_intensity: the intensity relative to the highest feature in the sample.
        area: the area of the peak
        rel_area: the area relative to the feature with the highest area in the sample.
        Spectrum: a matchms Spectrum object instance using data from msms
        samples: samples in which feature was detected.
        area_per_sample: list of SampleInfo instances summarizing area per sample
        height_per_sample: list of SampleInfo instances summarizing height per sample
        blank: bool to indicate if feature is blank-associated (if provided).
        groups: association to categories and groups is such data was provided.
        group_factors: indicates the group factors(fold differences) if provided.
        Annotations: objects summarizing associated annotation data
        networks: dict of objects representing associated networking data
        Scores: Object representing feature-associated scores
    """

    f_id: Optional[int] = None
    mz: Optional[float] = None
    rt: Optional[float] = None
    rt_start: Optional[float] = None
    rt_stop: Optional[float] = None
    rt_range: Optional[float] = None
    trace_rt: Optional[tuple] = None
    trace_int: Optional[tuple] = None
    fwhm: Optional[float] = None
    intensity: Optional[float] = None
    rel_intensity: Optional[float] = None
    area: Optional[float] = None
    rel_area: Optional[float] = None
    Spectrum: Optional[Any] = None
    samples: Optional[set] = None
    area_per_sample: Optional[list] = None
    height_per_sample: Optional[list] = None
    blank: Optional[bool] = None
    groups: Optional[dict] = None
    group_factors: Optional[dict] = None
    Annotations: Optional[Any] = None
    networks: Optional[dict] = None
    Scores: Optional[Any] = None

    def to_json(self: Self) -> dict:
        """Convert class attributes to json-compatible dict.

        Returns:
            A dictionary with class attributes as keys
        """
        attributes = (
            ("f_id", self.f_id, int),
            ("mz", self.mz, float),
            ("rt", self.rt, float),
            ("rt_start", self.rt_start, float),
            ("rt_stop", self.rt_stop, float),
            ("rt_range", self.rt_range, float),
            ("trace_rt", self.trace_rt, list),
            ("trace_int", self.trace_int, list),
            ("fwhm", self.fwhm, float),
            ("intensity", self.intensity, float),
            ("rel_intensity", self.rel_intensity, float),
            ("area", self.area, float),
            ("rel_area", self.rel_area, float),
            ("samples", self.samples, list),
            ("blank", self.blank, bool),
        )

        json_dict = {}
        for attribute in attributes:
            if attribute[1] is not None:
                json_dict[attribute[0]] = attribute[2](attribute[1])

        def _add_per_sample(attr: str):
            if getattr(self, attr) is not None:
                json_dict[attr] = [val.to_json() for val in getattr(self, attr)]
            else:
                json_dict[attr] = []

        _add_per_sample("area_per_sample")
        _add_per_sample("height_per_sample")

        if self.groups is not None and len(self.groups) != 0:
            json_dict["groups"] = {key: list(val) for key, val in self.groups.items()}
        else:
            json_dict["groups"] = {}

        if self.group_factors is not None:
            json_dict["group_factors"] = {
                key: [item.to_json() for item in val]
                for key, val in self.group_factors.items()
            }
        else:
            json_dict["group_factors"] = {}

        if self.Scores is not None:
            json_dict["scores"] = self.Scores.to_json()
        else:
            json_dict["scores"] = {}

        if self.Annotations is not None:
            json_dict["annotations"] = self.Annotations.to_json()
        else:
            json_dict["annotations"] = {}

        if self.networks is not None:
            json_dict["networks"] = {
                key: val.to_json() for key, val in self.networks.items()
            }
        else:
            json_dict["networks"] = {}

        json_dict["spectrum"] = {}

        return json_dict
