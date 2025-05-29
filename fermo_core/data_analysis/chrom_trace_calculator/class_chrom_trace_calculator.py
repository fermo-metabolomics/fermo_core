"""Calculate time/intensity traces of pseudo-chromatograms for plotting in GUI.

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
from typing import Optional, Self

from fermo_core.data_processing.builder_sample.dataclass_sample import Sample
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats

logger = logging.getLogger("fermo_core")


class ChromTraceData:
    """Internal class to organize attributes to create a pseudo-chromatogram trace

    Attributes:
        feature_id: id of molecular feature associated to data
        sample_id: id of sample associated to data
        rt_begin: RT of beginning of chromatogram
        int_begin: rel int of beginning of chromatogram
        rt_left_kink: RT of artificial point between start and fwhm for vis purpose
        int_left_kink: rel int of artificial point between start and fwhm for vis purpose
        rt_left_fwhm: RT of left side of feature width at half maximum
        int_left_fwhm: rel int of left side of feature width at half maximum
        rt_apex: RT of apex of the peak
        int_apex: rel int of apex of the peak
        rt_right_fwhm: rRT of right side of feature width at half maximum
        int_right_fwhm: rel int of right side of feature width at half max
        rt_right_kink: RT of artificial point between fwhm for vis purpose
        int_right_kink: rel int of artificial point between fwhm for vis purpose
        rt_end: RT of end of chromatogram
        int_end: rel int of end of chromatogram
        fwhm: feature width at half maximum intensity in minutes
        rt_range: peak range (from start to stop) in minutes
    """

    def __init__(
        self: Self,
        feature_id: int,
        sample_id: str,
        rt: float,
        rt_begin: float,
        rt_end: float,
    ):
        self.feature_id: int = feature_id
        self.sample_id: str = sample_id
        self.rt_begin: float = rt_begin
        self.int_begin: Optional[float] = None
        self.rt_left_kink: Optional[float] = None
        self.int_left_kink: Optional[float] = None
        self.rt_left_fwhm: Optional[float] = None
        self.int_left_fwhm: Optional[float] = None
        self.rt_apex: float = rt
        self.int_apex: Optional[float] = None
        self.rt_right_fwhm: Optional[float] = None
        self.int_right_fwhm: Optional[float] = None
        self.rt_right_kink: Optional[float] = None
        self.int_right_kink: Optional[float] = None
        self.rt_end: float = rt_end
        self.int_end: Optional[float] = None
        self.fwhm: Optional[float] = None
        self.rt_range: Optional[float] = None


class ChromTraceCalculator:
    """Organize methods to create a pseudo-chromatogram trace"""

    def __init__(self: Self):
        self.chrom_trace: Optional[ChromTraceData] = None

    def modify_samples(self: Self, samples: Repository, stats: Stats) -> Repository:
        """Load sample from repository for modification, then store again

        Attributes:
            samples: a Repository object containing the samples.
            stats: a Stats object containing a tuple of sample IDs.

        Returns:
            A Repository containing the modified samples.
        """
        logger.info(
            "'ChromTraceCalculator': started calculating pseudo-chromatogram traces."
        )

        for sample_id in stats.samples:
            sample = samples.get(sample_id)
            sample = self.modify_features_in_sample(sample)
            samples.modify(sample_id, sample)

        logger.info(
            "'ChromTraceCalculator': completed calculating pseudo-chromatogram traces."
        )
        return samples

    def modify_features_in_sample(self: Self, sample: Sample) -> Sample:
        """Add chromatogram traces to Features in Sample

        Attributes:
            sample: a Sample object containing Feature objects

        Returns:
            A Sample object with Feature objects with added chromatogram traces.
        """
        for feature in sample.feature_ids:
            self.chrom_trace = ChromTraceData(
                feature_id=feature,
                sample_id=sample.s_id,
                rt=sample.features.get(feature).rt,
                rt_begin=sample.features.get(feature).rt_start,
                rt_end=sample.features.get(feature).rt_stop,
            )
            self.validate_fwhm_rt(
                sample.features.get(feature).fwhm,
                sample.features.get(feature).rt_range,
            )
            self.calc_rt_left_fwhm()
            self.calc_rt_right_fwhm()
            self.calc_rt_left_kink()
            self.calc_rt_right_kink()
            self.assign_relative_intensity(
                sample.features.get(feature).rel_intensity,
            )

            sample.features.get(feature).trace_rt = self.create_trace_rt()
            sample.features.get(feature).trace_int = self.create_trace_int()

        return sample

    def validate_fwhm_rt(self: Self, fwhm: float, rt_range: float):
        """Validate that fwhm is less than rt_range and assign value.

        Attributes:
            fwhm: the full width at half maximum intensity in minutes
            rt_range: indicates range from beginning to end of peak, in minutes
        """
        if fwhm < rt_range:
            self.chrom_trace.fwhm = fwhm
            self.chrom_trace.rt_range = rt_range
        else:
            self.chrom_trace.fwhm = rt_range
            self.chrom_trace.rt_range = rt_range
            logger.debug(
                f"'ChromTraceCalculator': feature '{self.chrom_trace.feature_id}' in "
                f"sample "
                f"'{self.chrom_trace.sample_id}': "
                f"'fwhm' '{fwhm}' is >= 'rt_range' '{rt_range}'. "
                f"Setting 'fwhm' to 'rt_range'."
            )

    def calc_rt_left_fwhm(self: Self):
        """Calculate rt_left_fwhm, validate if in rt range"""
        rt_left_fwhm = self.chrom_trace.rt_apex - (0.5 * self.chrom_trace.fwhm)
        if rt_left_fwhm >= self.chrom_trace.rt_begin:
            self.chrom_trace.rt_left_fwhm = rt_left_fwhm
        else:
            self.chrom_trace.rt_left_fwhm = self.chrom_trace.rt_begin
            logger.debug(
                f"'ChromTraceCalculator': feature '{self.chrom_trace.feature_id}' in "
                f"sample "
                f"'{self.chrom_trace.sample_id}': "
                f"'rt_left_fwhm' '{rt_left_fwhm}' is < 'rt_begin' "
                f"'{self.chrom_trace.rt_begin}'. "
                f"Setting 'rt_left_fwhm' to 'rt_begin'."
            )

    def calc_rt_right_fwhm(self: Self):
        """Calculate rt_right_fwhm, validate if in rt range"""
        rt_right_fwhm = self.chrom_trace.rt_left_fwhm + self.chrom_trace.fwhm
        if rt_right_fwhm <= self.chrom_trace.rt_end:
            self.chrom_trace.rt_right_fwhm = rt_right_fwhm
        else:
            self.chrom_trace.rt_right_fwhm = self.chrom_trace.rt_end
            logger.debug(
                f"'ChromTraceCalculator': feature '{self.chrom_trace.feature_id}' in "
                f"sample "
                f"'{self.chrom_trace.sample_id}': "
                f"'rt_right_fwhm' '{rt_right_fwhm}' is > 'rt_end' "
                f"'{self.chrom_trace.rt_end}'. "
                f"Setting 'rt_right_fwhm' to 'rt_end'."
            )

    def calc_rt_left_kink(self: Self):
        """Calculate rt_left_kink and validate if in rt_range"""
        self.chrom_trace.rt_left_kink = self.chrom_trace.rt_left_fwhm - (
            0.5 * (self.chrom_trace.rt_left_fwhm - self.chrom_trace.rt_begin)
        )

    def calc_rt_right_kink(self: Self):
        """Calculate rt_right_kink and validate if in rt_range"""
        self.chrom_trace.rt_right_kink = self.chrom_trace.rt_right_fwhm + (
            0.5 * (self.chrom_trace.rt_end - self.chrom_trace.rt_right_fwhm)
        )

    def assign_relative_intensity(self: Self, relative_int: float):
        """Calculates and assigns relative intensity values

        Arguments:
            relative_int: intensity relative to the highest/most intense mol feature
        """
        self.chrom_trace.int_begin = relative_int * 0.0
        self.chrom_trace.int_left_kink = relative_int * 0.15
        self.chrom_trace.int_left_fwhm = relative_int * 0.5
        self.chrom_trace.int_apex = relative_int * 1.0
        self.chrom_trace.int_right_fwhm = relative_int * 0.5
        self.chrom_trace.int_right_kink = relative_int * 0.15
        self.chrom_trace.int_end = relative_int * 0.0

    def create_trace_rt(self: Self) -> tuple[float, ...]:
        """Create a trace of retention time points to assign to Feature object.

        Returns:
            A tuple of floats for RT of pseudo-chromatogram.
        """
        return (
            round(self.chrom_trace.rt_begin, 3),
            round(self.chrom_trace.rt_left_kink, 3),
            round(self.chrom_trace.rt_left_fwhm, 3),
            round(self.chrom_trace.rt_apex, 3),
            round(self.chrom_trace.rt_right_fwhm, 3),
            round(self.chrom_trace.rt_right_kink, 3),
            round(self.chrom_trace.rt_end, 3),
        )

    def create_trace_int(self: Self) -> tuple[float, ...]:
        """Create a trace of relative intensity points to assign to Feature object.

        Returns:
            A tuple of floats for int of pseudo-chromatogram.
        """
        return (
            round(self.chrom_trace.int_begin, 3),
            round(self.chrom_trace.int_left_kink, 3),
            round(self.chrom_trace.int_left_fwhm, 3),
            round(self.chrom_trace.int_apex, 3),
            round(self.chrom_trace.int_right_fwhm, 3),
            round(self.chrom_trace.int_right_kink, 3),
            round(self.chrom_trace.int_end, 3),
        )
