"""Organize the calling of annotation modules.

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
from pyexpat import features
from typing import Self

from pydantic import BaseModel

from fermo_core.config.class_default_settings import DefaultPaths
from fermo_core.data_analysis.annotation_manager.class_adduct_annotator import (
    AdductAnnotator,
)
from fermo_core.data_analysis.annotation_manager.class_fragment_annotator import (
    FragmentAnnotator,
)
from fermo_core.data_analysis.annotation_manager.class_mod_cos_annotator import (
    ModCosAnnotator,
)
from fermo_core.data_analysis.annotation_manager.class_ms2deepscore_annotator import (
    Ms2deepscoreAnnotator,
)
from fermo_core.data_analysis.annotation_manager.class_ms2query_annotator import (
    MS2QueryAnnotator,
)
from fermo_core.data_analysis.annotation_manager.class_mzmine_ann_parser import (
    MzmineAnnParser,
)
from fermo_core.data_analysis.annotation_manager.class_neutral_loss_annotator import (
    NeutralLossAnnotator,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager
from fermo_core.utils.utility_method_manager import UtilityMethodManager

logger = logging.getLogger("fermo_core")


class AnnotationManager(BaseModel):
    """Pydantic-based class to organize calling and logging of annotation modules

    Attributes:
        params: ParameterManager object, holds user-provided parameters
        stats: Stats object, holds stats on molecular features and samples
        features: Repository object, holds "General Feature" objects
        samples: Repository object, holds "Sample" objects
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository

    def return_attrs(
        self: Self,
    ) -> tuple[Stats, Repository, Repository, ParameterManager]:
        """Returns modified attributes from AnnotationManager to the calling function

        Returns:
            Tuple containing Stats, Feature Repository and Sample Repository objects.
        """
        return self.stats, self.features, self.samples, self.params

    def run_analysis(self: Self):
        """Organizes calling of data analysis steps."""
        logger.info("'AnnotationManager': started analysis steps.")

        def _eval_mzmine_file() -> bool:
            return (
                True
                if self.params.PeaktableParameters.format in ["mzmine3", "mzmine4"]
                else False
            )

        def _eval_ms2query_results_file() -> bool:
            return True if self.params.MS2QueryResultsParameters else False

        def _eval_as_results_file() -> bool:
            return True if self.params.AsResultsParameters else False

        modules = (
            (
                self.params.SpectralLibMatchingCosineParameters,
                self.run_user_lib_mod_cosine_matching,
            ),
            (
                self.params.SpectralLibMatchingDeepscoreParameters,
                self.run_user_lib_ms2deepscore_matching,
            ),
            (
                self.params.AdductAnnotationParameters,
                self.run_feature_adduct_annotation,
            ),
            (
                self.params.NeutralLossParameters,
                self.run_neutral_loss_annotation,
            ),
            (
                self.params.FragmentAnnParameters,
                self.run_fragment_annotation,
            ),
        )

        for module in modules:
            if getattr(module[0], "activate_module", False):
                module[1]()

        modules = (
            (
                _eval_mzmine_file(),
                self.run_mzmine_assignment,
            ),
            (
                _eval_ms2query_results_file(),
                self.run_ms2query_results_assignment,
            ),
            (
                _eval_as_results_file(),
                self.run_as_kcb_cosine_annotation,
            ),
            (
                _eval_as_results_file(),
                self.run_as_kcb_deepscore_annotation,
            ),
        )
        for module in modules:
            if module[0]:
                module[1]()

        logger.info("'AnnotationManager': completed analysis steps.")

    def verify_user_lib_params(self: Self) -> bool:
        """Perform preliminary checks on user-provided library and params.

        Returns:
            A bool indicating pass or fail
        """
        if self.params.SpecLibParameters is None:
            logger.warning(
                "'AnnotationManager': no spectral library params provided - SKIP"
            )
            return False
        elif self.stats.spectral_library is None:
            logger.warning(
                "'AnnotationManager': no spectral library file provided - SKIP"
            )
            return False
        elif len(self.stats.spectral_library) == 0:
            logger.warning("'AnnotationManager': spectral library file is empty - SKIP")
            return False
        else:
            return True

    def run_user_lib_mod_cosine_matching(self: Self):
        """Match features against a user-provided spectral library using mod cosine."""
        logger.info(
            "'AnnotationManager': started matching of features against a "
            "user-provided spectral library using the modified cosine algorithm."
        )

        if not self.verify_user_lib_params():
            return

        try:
            mod_cosine_annotator = ModCosAnnotator(
                features=self.features,
                active_features=self.stats.active_features,
                library=self.stats.spectral_library,
                library_name=self.params.SpecLibParameters.dirpath.name,
                fragment_tol=self.params.SpectralLibMatchingCosineParameters.fragment_tol,
                score_cutoff=self.params.SpectralLibMatchingCosineParameters.score_cutoff,
                min_nr_matched_peaks=self.params.SpectralLibMatchingCosineParameters.min_nr_matched_peaks,
                max_precursor_mass_diff=self.params.SpectralLibMatchingCosineParameters.max_precursor_mass_diff,
            )
            mod_cosine_annotator.prepare_queries()
            mod_cosine_annotator.calculate_scores_mod_cosine()
            mod_cosine_annotator.extract_userlib_scores()
            self.features = mod_cosine_annotator.return_features()
            self.params.SpectralLibMatchingCosineParameters.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of ModCosineAnnotator - SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed matching of features against a "
            "user-provided spectral library using the modified cosine algorithm."
        )

    def run_user_lib_ms2deepscore_matching(self: Self):
        """Match features against user-provided spectral library using ms2deepscore."""
        logger.info(
            "'AnnotationManager': started matching of features against a "
            "user-provided spectral library using the ms2deepscore algorithm."
        )

        if not self.verify_user_lib_params():
            return

        try:
            ms2deepscore_annotator = Ms2deepscoreAnnotator(
                features=self.features,
                active_features=self.stats.active_features,
                polarity=self.params.PeaktableParameters.polarity,
                library=self.stats.spectral_library,
                library_name=self.params.SpecLibParameters.dirpath.name,
                score_cutoff=self.params.SpectralLibMatchingDeepscoreParameters.score_cutoff,
                max_precursor_mass_diff=self.params.SpectralLibMatchingDeepscoreParameters.max_precursor_mass_diff,
            )
            ms2deepscore_annotator.prepare_queries()
            ms2deepscore_annotator.calculate_scores_ms2deepscore()
            ms2deepscore_annotator.extract_userlib_scores()
            self.features = ms2deepscore_annotator.return_features()
            self.params.SpectralLibMatchingDeepscoreParameters.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of Ms2deepscoreAnnotator "
                "- SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed matching of features against a "
            "user-provided spectral library using the ms2deepscore algorithm."
        )

    def run_feature_adduct_annotation(self: Self):
        """Perform feature adduct annotation"""
        logger.info("'AnnotationManager': started feature adduct annotation.")

        try:
            adduct_annotator = AdductAnnotator(
                params=self.params,
                stats=self.stats,
                features=self.features,
                samples=self.samples,
            )
            adduct_annotator.run_analysis()
            self.features, self.stats = adduct_annotator.return_features()
            self.params.AdductAnnotationParameters.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of AdductAnnotator " "- SKIP"
            )
            return

        logger.info("'AnnotationManager': completed feature adduct annotation.")

    def run_neutral_loss_annotation(self: Self):
        """Perform feature MS2 neutral loss annotation"""
        logger.info("'AnnotationManager': started feature neutral loss annotation.")

        try:
            neutralloss_annotator = NeutralLossAnnotator(
                params=self.params,
                stats=self.stats,
                features=self.features,
                samples=self.samples,
            )
            neutralloss_annotator.run_analysis()
            self.features, self.params = neutralloss_annotator.return_attributes()
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of NeutralLossAnnotator "
                "- SKIP"
            )
            return

        logger.info("'AnnotationManager': completed feature neutral loss annotation.")

    def run_fragment_annotation(self: Self):
        """Perform feature MS2 fragment annotation"""
        logger.info("'AnnotationManager': started feature fragment annotation.")

        try:
            fragment_annotator = FragmentAnnotator(
                params=self.params,
                stats=self.stats,
                features=self.features,
                samples=self.samples,
            )
            fragment_annotator.run_analysis()
            self.features, self.params = fragment_annotator.return_attributes()
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of FragmentAnnotator "
                "- SKIP"
            )
            return

        logger.info("'AnnotationManager': completed feature fragment annotation.")

    def run_mzmine_assignment(self) -> None:
        """Annotate Features from optional mzmine file fields"""
        logger.info(
            "'AnnotationManager': started annotation from existing Mzmine peaktable."
        )

        try:
            mzmine_ann = MzmineAnnParser(
                params=self.params, stats=self.stats, features=self.features
            )
            mzmine_ann.run()
            self.features = mzmine_ann.return_attributes()
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during Mzmine peaktable annotation assignment - SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed annotation from existing Mzmine peaktable."
        )

    def run_ms2query_results_assignment(self: Self):
        """Annotate Features from existing MS2Query results"""

        logger.info(
            "'AnnotationManager': started annotation from existing MS2Query results."
        )

        try:
            ms2query_annotator = MS2QueryAnnotator(
                params=self.params,
                features=self.features,
                active_features=self.stats.active_features,
                cutoff=self.params.MS2QueryResultsParameters.score_cutoff,
            )
            ms2query_annotator.assign_feature_info(
                self.params.MS2QueryResultsParameters.filepath,
            )
            self.features = ms2query_annotator.return_features()

        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during MS2Query Results Assignment - SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed annotation from existing MS2Query results."
        )

    def run_as_kcb_cosine_annotation(self: Self):
        """Match features against a antiSMASH knownclusterblast-derived library.

        Allows modified cosine-based library matching against the in silico generated
        MS2 spectra of significant antiSMASH KnownClusterBlast (MIBiG) matches.
        """
        logger.info(
            "'AnnotationManager': started antiSMASH KnownClusterBlast "
            "modified cosine annotation."
        )

        if self.params.PeaktableParameters.polarity == "negative":
            logger.warning(
                "'AnnotationManager': negative ion mode detected. antiSMASH "
                "KnownClusterBlast result annotation only available for positive ion "
                "mode - SKIP"
            )
            return

        if not (
            self.params.AsKcbCosineMatchingParams
            and self.params.AsKcbCosineMatchingParams.activate_module
        ):
            logger.warning(
                "'AnnotationManager': antiSMASH results file provided but "
                "'as_kcb_matching/modified_cosine' is turned off - SKIP"
            )
            return

        try:
            kcb_results = UtilityMethodManager().extract_as_kcb_results(
                as_results=self.params.AsResultsParameters.directory_path,
                cutoff=self.params.AsResultsParameters.similarity_cutoff,
            )
            mibig_bgcs = {key for key, value in kcb_results.items()}
            spec_library = UtilityMethodManager().create_mibig_spec_lib(mibig_bgcs)
            kcb_annotator = ModCosAnnotator(
                features=self.features,
                active_features=self.stats.active_features,
                library=spec_library,
                library_name=DefaultPaths().library_mibig_pos.name,
                fragment_tol=self.params.AsKcbCosineMatchingParams.fragment_tol,
                score_cutoff=self.params.AsKcbCosineMatchingParams.score_cutoff,
                min_nr_matched_peaks=self.params.AsKcbCosineMatchingParams.min_nr_matched_peaks,
                max_precursor_mass_diff=self.params.AsKcbCosineMatchingParams.max_precursor_mass_diff,
            )
            kcb_annotator.prepare_queries()
            kcb_annotator.calculate_scores_mod_cosine()
            kcb_annotator.extract_mibig_scores(kcb_results)
            self.features = kcb_annotator.return_features()
            self.params.AsKcbCosineMatchingParams.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of antiSMASH "
                "KnownClusterBlast modified cosine annotation - SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed antiSMASH KnownClusterBlast "
            "modified cosine annotation."
        )

    def run_as_kcb_deepscore_annotation(self: Self):
        """Match features against a antiSMASH knownclusterblast-derived library.

        Allows MS2DeepScore-based library matching against the in silico generated
        MS2 spectra of significant antiSMASH KnownClusterBlast (MIBiG) matches.
        """
        logger.info(
            "'AnnotationManager': started antiSMASH KnownClusterBlast MS2DeepScore "
            "annotation."
        )

        if self.params.PeaktableParameters.polarity == "negative":
            logger.warning(
                "'AnnotationManager': negative ion mode detected. antiSMASH "
                "KnownClusterBlast result annotation only available for positive ion "
                "mode - SKIP"
            )
            return

        if not (
            self.params.AsKcbDeepscoreMatchingParams
            and self.params.AsKcbDeepscoreMatchingParams.activate_module
        ):
            logger.warning(
                "'AnnotationManager': antiSMASH results file provided but "
                "'as_kcb_matching/ms2deepscore' is turned off - SKIP"
            )
            return

        try:
            kcb_results = UtilityMethodManager().extract_as_kcb_results(
                as_results=self.params.AsResultsParameters.directory_path,
                cutoff=self.params.AsResultsParameters.similarity_cutoff,
            )
            mibig_bgcs = {key for key, value in kcb_results.items()}
            spec_library = UtilityMethodManager().create_mibig_spec_lib(mibig_bgcs)

            kcb_annotator = Ms2deepscoreAnnotator(
                features=self.features,
                active_features=self.stats.active_features,
                polarity=self.params.PeaktableParameters.polarity,
                library=spec_library,
                library_name=DefaultPaths().library_mibig_pos.name,
                score_cutoff=self.params.AsKcbDeepscoreMatchingParams.score_cutoff,
                max_precursor_mass_diff=self.params.AsKcbDeepscoreMatchingParams.max_precursor_mass_diff,
            )
            kcb_annotator.prepare_queries()
            kcb_annotator.calculate_scores_ms2deepscore()
            kcb_annotator.extract_mibig_scores(kcb_results)
            self.features = kcb_annotator.return_features()
            self.params.AsKcbDeepscoreMatchingParams.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'AnnotationManager': Error during running of antiSMASH "
                "KnownClusterBlast MS2DeepScore annotation - SKIP"
            )
            return

        logger.info(
            "'AnnotationManager': completed antiSMASH KnownClusterBlast MS2DeepScore "
            "annotation."
        )
