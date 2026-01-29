"""Organizes classes that hold and validate parameters.

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

from pydantic import (
    BaseModel,
    DirectoryPath,
    FilePath,
    PositiveFloat,
    PositiveInt,
    model_validator,
)

from fermo_core.input_output.class_validation_manager import ValidationManager

logger = logging.getLogger("fermo_core")


class PeaktableParameters(BaseModel):
    """A Pydantic-based class for representing and validating peaktable parameters.

    Attributes:
        filepath: a pathlib Path object pointing towards a peaktable file
        format: indicates the format of the peaktable file
        polarity: indicates the polarity of the data ('positive', 'negative').

    Raise:
        ValueError: Unsupported peaktable format.
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    filepath: FilePath
    format: str
    polarity: str

    @model_validator(mode="after")
    def val(self):
        if self.format == "mzmine3" or self.format == "mzmine4":
            ValidationManager.validate_file_extension(self.filepath, ".csv")
            ValidationManager.validate_csv_file(self.filepath)
            ValidationManager.validate_csv_has_rows(self.filepath)
            ValidationManager.validate_peaktable_mzmine(self.filepath)
            ValidationManager.validate_no_duplicate_entries_csv_column(
                self.filepath, "id"
            )
        else:
            raise ValueError(f"Unsupported peaktable format: '{self.format}'.")
        ValidationManager.validate_allowed(self.polarity, ["positive", "negative"])
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "filepath": str(self.filepath.name),
            "format": str(self.format),
            "polarity": str(self.polarity),
        }


class MsmsParameters(BaseModel):
    """A Pydantic-based class for representing and validating MS/MS file parameters.

    Attributes:
        filepath: a pathlib Path object pointing towards an MS/MS file
        format: indicates the format of the MS/MS file
        rel_int_from: the minimum relative intensity of MS2 fragments to be retained

    Raise:
        ValueError: Unsupported MS/MS file format.
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    filepath: FilePath
    format: str
    rel_int_from: float

    @model_validator(mode="after")
    def val(self):
        match self.format:
            case "mgf":
                ValidationManager.validate_file_extension(self.filepath, ".mgf")
                ValidationManager.validate_mgf_file(self.filepath)
            case _:
                raise ValueError(f"Unsupported MS/MS format: '{self.format}'.")
        ValidationManager.validate_float_zero_one(self.rel_int_from)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "filepath": str(self.filepath.name),
            "format": self.format,
            "rel_int_from": self.rel_int_from,
        }


class PhenotypeParameters(BaseModel):
    """A Pydantic-based class for representing and validating phenotype file parameters.

    Attributes:
        filepath: a pathlib Path object pointing towards a phenotype file
        format: indicates the format of the phenotype file

    Raise:
        ValueError: Unsupported phenotype file format.
    """

    filepath: FilePath
    format: str

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_file_extension(self.filepath, ".csv")
        ValidationManager.validate_csv_file(self.filepath)
        ValidationManager.validate_csv_has_rows(self.filepath)
        match self.format:
            case "qualitative":
                ValidationManager.validate_pheno_qualitative(self.filepath)
                ValidationManager.validate_no_duplicate_entries_csv_column(
                    self.filepath, "sample_name"
                )
            case "quantitative-percentage":
                ValidationManager.validate_pheno_quant_percentage(self.filepath)
                ValidationManager.validate_no_duplicate_entries_csv_column(
                    self.filepath, "well"
                )
            case "quantitative-concentration":
                ValidationManager.validate_pheno_quant_concentration(self.filepath)
                ValidationManager.validate_no_duplicate_entries_csv_column(
                    self.filepath, "well"
                )
            case _:
                raise ValueError(f"Unsupported phenotype format: '{self.format}'.")
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "filepath": str(self.filepath.name),
            "format": str(self.format),
        }


class GroupMetadataParameters(BaseModel):
    """A Pydantic-based class for representing and validating group metadata parameters.

    Attributes:
        filepath: a pathlib Path object pointing towards a group metadata file
        format: indicates the format of the group metadata file

    Raise:
        ValueError: Unsupported group metadata file format.
        pydantic.ValidationError: Pydantic validation failed during instantiation.

    Notes:
        A fermo-style group data file is a .csv-file with the layout:

        sample_name,group_col_1,group_col_2,...,group_col_n \n
        sample1,medium_A,condition_A \n
        sample2,medium_B,condition_A\n
        sample3,medium_C,condition_A \n

        Ad values: The only prohibited value is 'DEFAULT' which is reserved for
        internal use. 'BLANK' os a special value that indicates the sample/medium
        blank for automated subtraction.
    """

    filepath: FilePath
    format: str

    @model_validator(mode="after")
    def val(self):
        match self.format:
            case "fermo":
                ValidationManager.validate_file_extension(self.filepath, ".csv")
                ValidationManager.validate_csv_file(self.filepath)
                ValidationManager.validate_csv_has_rows(self.filepath)
                ValidationManager.validate_group_metadata_fermo(self.filepath)
                ValidationManager.validate_no_duplicate_entries_csv_column(
                    self.filepath, "sample_name"
                )
            case _:
                raise ValueError(
                    f"Unsupported group metadata format: " f"'{self.format}'."
                )
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "filepath": str(self.filepath.name),
            "format": str(self.format),
        }


class SpecLibParameters(BaseModel):
    """Pydantic-based class for repres. and valid. of spectral library files parameters.

    Attributes:
        dirpath: a pathlib Path object pointing towards a dir containing spec lib files
        format: indicates the format of the spectral library files

    Raise:
        ValueError: Unsupported spectral library format.
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    dirpath: DirectoryPath
    format: str

    @model_validator(mode="after")
    def val(self):
        match self.format:
            case "mgf":
                for f in self.dirpath.iterdir():
                    ValidationManager.validate_file_extension(f, ".mgf")
                    ValidationManager.validate_mgf_file(f)
            case _:
                raise ValueError(f"Unsupported library format '{self.format}'.")
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "dirpath": str(self.dirpath.name),
            "format": str(self.format),
        }


class MS2QueryResultsParameters(BaseModel):
    """Pydantic-based class for repres. and valid. of MS2Query result parameters.

    Attributes:
        filepath: a pathlib Path object pointing towards a MS2Query results file
        score_cutoff: the minimal score to retain the annotation

    Raise:
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    filepath: FilePath
    score_cutoff: PositiveFloat

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_file_extension(self.filepath, ".csv")
        ValidationManager.validate_csv_has_rows(self.filepath)
        ValidationManager.validate_ms2query_results(self.filepath)
        ValidationManager.validate_float_zero_one(self.score_cutoff)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "filepath": str(self.filepath.name),
            "score_cutoff": self.score_cutoff,
        }


class AsResultsParameters(BaseModel):
    """A Pydantic-based class for representing and validating an antiSMASH results dir.

    Attributes:
        directory_path: the output directory path
        similarity_cutoff: a fraction indicating the minimum shared similarity required

    Raise:
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    directory_path: DirectoryPath
    similarity_cutoff: PositiveFloat

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_float_zero_one(self.similarity_cutoff)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        return {
            "directory_path": str(self.directory_path.name),
            "similarity_cutoff": self.similarity_cutoff,
        }


class OutputParameters(BaseModel):
    """A Pydantic-based class for representing and validating output parameters.

    Attributes:
        directory_path: the output directory path
    """

    directory_path: DirectoryPath

    @model_validator(mode="after")
    def val(self):
        self.directory_path = self.directory_path.joinpath("results")
        self.directory_path.mkdir(exist_ok=True)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        try:
            return {
                "directory_path": f"{self.directory_path.parent.name}/{self.directory_path.name}"
            }
        except AttributeError:
            return {"directory_path": "not specified"}


class AdductAnnotationParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of adduct annotation parameters.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        mass_dev_ppm: The estimated maximum mass deviation in ppm.
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    mass_dev_ppm: PositiveFloat
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_mass_deviation_ppm(self.mass_dev_ppm)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "mass_dev_ppm": float(self.mass_dev_ppm),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class NeutralLossParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of neutral loss annotation params.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        mass_dev_ppm: The estimated maximum mass deviation in ppm.
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    mass_dev_ppm: PositiveFloat
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_mass_deviation_ppm(self.mass_dev_ppm)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "mass_dev_ppm": float(self.mass_dev_ppm),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class FragmentAnnParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of fragment annotation params.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        mass_dev_ppm: The estimated maximum mass deviation in ppm.
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    mass_dev_ppm: PositiveFloat
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_mass_deviation_ppm(self.mass_dev_ppm)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "mass_dev_ppm": float(self.mass_dev_ppm),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class SpecSimNetworkCosineParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of spec. similarity network params.

    This class manages parameters for the modified cosine algorithm.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        msms_min_frag_nr: minimum nr of fragments a spectrum must have to be considered.
        fragment_tol: the tolerance between matched fragments, in m/z units.
        score_cutoff: the minimum similarity score between two spectra.
        max_nr_links: max nr of connections from a node.
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    msms_min_frag_nr: PositiveInt
    fragment_tol: PositiveFloat
    score_cutoff: PositiveFloat
    max_nr_links: PositiveInt
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_float_zero_one(self.score_cutoff)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "msms_min_frag_nr": int(self.msms_min_frag_nr),
                "fragment_tol": float(self.fragment_tol),
                "score_cutoff": float(self.score_cutoff),
                "max_nr_links": int(self.max_nr_links),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class SpecSimNetworkDeepscoreParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of spec. similarity network params.

    This class manages parameters for the ms2deepscore algorithm.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        score_cutoff: the minimum similarity score between two spectra.
        max_nr_links: max links to a single spectra.
        msms_min_frag_nr: minimum number of fragments in MS2 to run it in analysis
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = True
    score_cutoff: PositiveFloat
    max_nr_links: PositiveInt
    msms_min_frag_nr: PositiveInt
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_float_zero_one(self.score_cutoff)
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "score_cutoff": float(self.score_cutoff),
                "max_nr_links": int(self.max_nr_links),
                "msms_min_frag_nr": int(self.msms_min_frag_nr),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class FeatureFilteringParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of feature filtering parameters.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        filter_rel_int_range_min: min value to filter feature for rel int
        filter_rel_int_range_max: max value to filter feature for rel int
        filter_rel_area_range_min: min value to filter feature for rel area
        filter_rel_area_range_max: max value to filter feature for rel area
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    filter_rel_int_range_min: float
    filter_rel_int_range_max: float
    filter_rel_area_range_min: float
    filter_rel_area_range_max: float
    module_passed: bool = False

    @model_validator(mode="after")
    def validate_attrs(self):
        ValidationManager.validate_range_zero_one(
            [self.filter_rel_int_range_min, self.filter_rel_int_range_max]
        )
        ValidationManager.validate_range_zero_one(
            [self.filter_rel_area_range_min, self.filter_rel_area_range_max]
        )
        return self

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "filter_rel_int_range_min": self.filter_rel_int_range_min,
                "filter_rel_int_range_max": self.filter_rel_int_range_max,
                "filter_rel_area_range_min": self.filter_rel_area_range_min,
                "filter_rel_area_range_max": self.filter_rel_area_range_max,
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class BlankAssignmentParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of blank assignment parameters.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        factor: An integer fold-change to differentiate blank features.
        algorithm: the algorithm to summarize values of different samples.
        value: the type of value to use for determination
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    factor: PositiveInt
    algorithm: str
    value: str
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_allowed(
            self.algorithm, ["mean", "median", "maximum"]
        )
        ValidationManager.validate_allowed(self.value, ["height", "area"])
        return self

    def to_json(self: Self) -> dict:
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "factor": int(self.factor),
                "algorithm": str(self.algorithm),
                "value": str(self.value),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class GroupFactAssignmentParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of group factor assignment params.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        algorithm: the algorithm to summarize values of different samples.
        value: the type of value to use for comparison
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    algorithm: str
    value: str
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_allowed(
            self.algorithm, ["mean", "median", "maximum"]
        )
        ValidationManager.validate_allowed(self.value, ["height", "area"])
        return self

    def to_json(self: Self) -> dict:
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "algorithm": str(self.algorithm),
                "value": str(self.value),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class PhenoQualAssgnParams(BaseModel):
    """A Pydantic-based class for phenotype qualitative assignment parameters

    Attributes:
        activate_module: bool to indicate if module should be executed.
        factor: An integer fold-change to differentiate phenotype-assoc. features.
        algorithm: the algorithm to summarize values of active vs inactive samples.
        value: the type of value to use for determination
        p_val_cutoff: minimum significance threshold derived by test; 0 to disable test
        test: type of test performed
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    factor: PositiveInt
    algorithm: str
    value: str
    p_val_cutoff: float = 0.0
    test: str = "None"
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_allowed(self.algorithm, ["mean", "median", "minmax"])
        ValidationManager.validate_allowed(self.value, ["height", "area"])
        ValidationManager.validate_allowed(self.test, ["None", "Welsh", "Wilcoxon"])
        ValidationManager.validate_float_zero_one(self.p_val_cutoff)
        return self

    def to_json(self: Self) -> dict:
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "factor": int(self.factor),
                "algorithm": str(self.algorithm),
                "value": str(self.value),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class PhenoQuantPercentAssgnParams(BaseModel):
    """A Pydantic-based class for phenotype quantitative percentage assignment params

    Attributes:
        activate_module: bool to indicate if module should be executed.
        sample_avg: algorithm to summarize mult measurements per sample for single assay
        value: the type of value to use for determination
        algorithm: the statistical algorithm to calculate correlation
        fdr_corr: false-discovery rate correction algorithm (one of statsmodels.stats.multitest.multipletests)
        p_val_cutoff: minimum FDR-corrected p-value to consider in assignment
        coeff_cutoff: minimum correlation coefficient cutoff to consider in assignment
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    sample_avg: str
    value: str
    algorithm: str
    fdr_corr: str
    p_val_cutoff: PositiveFloat
    coeff_cutoff: PositiveFloat
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_allowed(self.sample_avg, ["mean", "median"])
        ValidationManager.validate_allowed(self.value, ["area"])
        ValidationManager.validate_allowed(
            self.algorithm, ["pearson", "spearman", "spearman_permutation"]
        )
        ValidationManager.validate_allowed(
            self.fdr_corr,
            [
                "bonferroni",
                "sidak",
                "holm-sidak",
                "holm",
                "simes-hochberg",
                "hommel",
                "fdr_bh",
                "fdr_by",
                "fdr_tsbh",
                "fdr_tsbky",
            ],
        )
        ValidationManager.validate_float_zero_one(self.p_val_cutoff)
        ValidationManager.validate_float_zero_one(self.coeff_cutoff)
        return self

    def to_json(self: Self) -> dict:
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "sample_avg": self.sample_avg,
                "value": self.value,
                "algorithm": self.algorithm,
                "p_val_cutoff": self.p_val_cutoff,
                "fdr_corr": self.fdr_corr,
                "coeff_cutoff": self.coeff_cutoff,
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class PhenoQuantConcAssgnParams(BaseModel):
    """A Pydantic-based class for phenotype quantitative concentration assignment params

    Attributes:
        activate_module: bool to indicate if module should be executed.
        sample_avg: algorithm to summarize mult measurements per sample for single assay
        value: the type of value to use for determination
        algorithm: the statistical algorithm to calculate correlation
        fdr_corr: false-discovery rate correction algorithm (one of statsmodels.stats.multitest.multipletests)
        p_val_cutoff: minimum FDR-corrected p-value to consider in assignment
        coeff_cutoff: minimum correlation coefficient cutoff to consider in assignment
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    sample_avg: str
    value: str
    algorithm: str
    fdr_corr: str
    p_val_cutoff: PositiveFloat
    coeff_cutoff: PositiveFloat
    module_passed: bool = False

    @model_validator(mode="after")
    def val(self):
        ValidationManager.validate_allowed(self.sample_avg, ["mean", "median"])
        ValidationManager.validate_allowed(self.value, ["area"])
        ValidationManager.validate_allowed(
            self.algorithm, ["pearson", "spearman", "spearman_permutation"]
        )
        ValidationManager.validate_allowed(
            self.fdr_corr,
            [
                "bonferroni",
                "sidak",
                "holm-sidak",
                "holm",
                "simes-hochberg",
                "hommel",
                "fdr_bh",
                "fdr_by",
                "fdr_tsbh",
                "fdr_tsbky",
            ],
        )
        ValidationManager.validate_float_zero_one(self.p_val_cutoff)
        ValidationManager.validate_float_zero_one(self.coeff_cutoff)
        return self

    def to_json(self: Self) -> dict:
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "sample_avg": self.sample_avg,
                "value": self.value,
                "algorithm": self.algorithm,
                "fdr_corr": self.fdr_corr,
                "p_val_cutoff": self.p_val_cutoff,
                "coeff_cutoff": self.coeff_cutoff,
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class SpectralLibMatchingCosineParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of spectral library matching params.

    This class addresses parameters for the modified cosine algorithm.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        fragment_tol: max tolerable diff to consider two fragments as equal, in m/z
        min_nr_matched_peaks: peak cutoff to consider a match of two MS/MS spectra
        score_cutoff: score cutoff to consider a match of two MS/MS spectra
        max_precursor_mass_diff: maximum precursor mass difference
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    fragment_tol: PositiveFloat
    min_nr_matched_peaks: PositiveInt
    score_cutoff: PositiveFloat
    max_precursor_mass_diff: PositiveInt
    module_passed: bool = False

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "fragment_tol": float(self.fragment_tol),
                "min_nr_matched_peaks": int(self.min_nr_matched_peaks),
                "score_cutoff": float(self.score_cutoff),
                "max_precursor_mass_diff": int(self.max_precursor_mass_diff),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class SpectralLibMatchingDeepscoreParameters(BaseModel):
    """A Pydantic-based class for repr. and valid. of spectral library matching params.

    This class addresses parameters for the MS2DeepScore algorithm.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        score_cutoff: score cutoff to consider a match of two MS/MS spectra.
        max_precursor_mass_diff: max allowed precursor mz difference to accept a match
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    score_cutoff: PositiveFloat
    max_precursor_mass_diff: PositiveInt
    module_passed: bool = False

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "score_cutoff": float(self.score_cutoff),
                "max_precursor_mass_diff": int(self.max_precursor_mass_diff),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class AsKcbCosineMatchingParams(BaseModel):
    """Pydantic-based class for params of antiSMASH KCB results mod cosine matching.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        fragment_tol: max tolerable diff to consider two fragments as equal, in m/z
        min_nr_matched_peaks: peak cutoff to consider a match of two MS/MS spectra
        score_cutoff: score cutoff to consider a match of two MS/MS spectra
        max_precursor_mass_diff: maximum precursor mass difference
        module_passed: indicates that the module ran without errors
    """

    activate_module: bool = False
    fragment_tol: PositiveFloat
    min_nr_matched_peaks: PositiveInt
    score_cutoff: PositiveFloat
    max_precursor_mass_diff: PositiveInt
    module_passed: bool = False

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "fragment_tol": float(self.fragment_tol),
                "min_nr_matched_peaks": int(self.min_nr_matched_peaks),
                "score_cutoff": float(self.score_cutoff),
                "max_precursor_mass_diff": int(self.max_precursor_mass_diff),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}


class AsKcbDeepscoreMatchingParams(BaseModel):
    """Pydantic-based class for params of antiSMASH KCB results deepscore matching.

    Attributes:
        activate_module: bool to indicate if module should be executed.
        score_cutoff: score cutoff to consider a match of two MS/MS spectra.
        max_precursor_mass_diff: max allowed precursor mz difference to accept a match
        module_passed: indicates that the module ran without errors

    Raise:
        pydantic.ValidationError: Pydantic validation failed during instantiation.
    """

    activate_module: bool = False
    score_cutoff: PositiveFloat
    max_precursor_mass_diff: PositiveInt
    module_passed: bool = False

    def to_json(self: Self) -> dict:
        """Convert attributes to json-compatible ones."""
        if self.activate_module:
            return {
                "activate_module": self.activate_module,
                "score_cutoff": float(self.score_cutoff),
                "max_precursor_mass_diff": int(self.max_precursor_mass_diff),
                "module_passed": self.module_passed,
            }
        else:
            return {"activate_module": self.activate_module}
