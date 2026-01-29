"""Organize the calling of phenotype annotation modules.

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

from pydantic import BaseModel

from fermo_core.data_analysis.phenotype_manager.class_phen_qual_assigner import (
    PhenQualAssigner,
)
from fermo_core.data_analysis.phenotype_manager.class_phen_quant_ass import PhenQuantAss
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager

logger = logging.getLogger("fermo_core")


class PhenotypeManager(BaseModel):
    """Pydantic-based class to organize calling and logging of phenotype annot modules

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

    def return_attrs(self: Self) -> tuple[Stats, Repository, ParameterManager]:
        """Returns modified attributes from PhenotypeManager to the calling function

        Returns:
            Tuple containing Stats, Feature Repository objects.
        """
        return self.stats, self.features, self.params

    def run_analysis(self: Self):
        """Organizes calling of phenotype annotation steps."""
        logger.info("'PhenotypeManager': started analysis steps.")

        match self.params.PhenotypeParameters.format:
            case "qualitative":
                self.run_assigner_qualitative()
            case "quantitative-percentage":
                self.run_assigner_quant_percentage()
            case "quantitative-concentration":
                self.run_assigner_quant_concentration()
            case _:
                logger.warning("'PhenotypeManager': unexpected phenotype format - SKIP")
                return

        logger.info("'PhenotypeManager': completed analysis steps.")

    def run_assigner_qualitative(self: Self):
        """Run the phenotype feature annotation based on qualitative data"""
        logger.info("'PhenotypeManager': started qualitative phenotype data analysis.")

        if self.params.PhenoQualAssgnParams.activate_module is False:
            logger.info(
                "'PhenotypeManager': parameters for "
                "'phenotype_assignment/qualitative' not specified or module turned off "
                " - SKIP"
            )
            return

        try:
            qual_assigner = PhenQualAssigner(
                params=self.params,
                features=self.features,
                stats=self.stats,
                samples=self.samples,
            )
            qual_assigner.run_analysis()
            self.stats, self.features = qual_assigner.return_values()
            self.params.PhenoQualAssgnParams.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'PhenotypeManager': Error during running of PhenQualAssigner - SKIP"
            )
            return

        logger.info(
            "'PhenotypeManager': completed qualitative phenotype data analysis."
        )

    def run_assigner_quant_percentage(self: Self):
        """Run the phenotype feature annotation based on quantitative percentage data"""
        logger.info(
            "'PhenotypeManager': started quantitative phenotype data analysis for "
            "percentage data."
        )

        if self.params.PhenoQuantPercentAssgnParams.activate_module is False:
            logger.info(
                "'PhenotypeManager': parameters for "
                "'phenotype_assignment/quantitative-percentage' not specified or "
                "module turned off - SKIP"
            )
            return

        try:
            quant_p_assigner = PhenQuantAss(
                coeff_cutoff=self.params.PhenoQuantPercentAssgnParams.coeff_cutoff,
                p_val_cutoff=self.params.PhenoQuantPercentAssgnParams.p_val_cutoff,
                fdr_corr=self.params.PhenoQuantPercentAssgnParams.fdr_corr,
                mode="percentage",
                alg=self.params.PhenoQuantPercentAssgnParams.algorithm,
                features=self.features,
                stats=self.stats,
                samples=self.samples,
            )
            quant_p_assigner.run_analysis()
            self.stats, self.features = quant_p_assigner.return_values()
            self.params.PhenoQuantPercentAssgnParams.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'PhenotypeManager': Error during running of PhenQuantPercAssigner - SKIP"
            )
            return

        logger.info(
            "'PhenotypeManager': completed quantitative phenotype data analysis for "
            "percentage data."
        )

    def run_assigner_quant_concentration(self: Self):
        """Run the phenotype feature annotation based on quantitative concentr. data"""
        logger.info(
            "'PhenotypeManager': started quantitative phenotype data analysis for "
            "concentration data."
        )

        if self.params.PhenoQuantConcAssgnParams.activate_module is False:
            logger.info(
                "'PhenotypeManager': parameters for "
                "'phenotype_assignment/quantitative-concentration' not specified or "
                "module turned off - SKIP"
            )
            return

        try:
            quant_conc_assigner = PhenQuantAss(
                coeff_cutoff=self.params.PhenoQuantConcAssgnParams.coeff_cutoff,
                p_val_cutoff=self.params.PhenoQuantConcAssgnParams.p_val_cutoff,
                fdr_corr=self.params.PhenoQuantConcAssgnParams.fdr_corr,
                alg=self.params.PhenoQuantConcAssgnParams.algorithm,
                mode="percentage",
                features=self.features,
                stats=self.stats,
                samples=self.samples,
            )
            quant_conc_assigner.run_analysis()
            self.stats, self.features = quant_conc_assigner.return_values()
            self.params.PhenoQuantConcAssgnParams.module_passed = True
        except Exception as e:
            logger.error(str(e))
            logger.error(
                "'PhenotypeManager': Error during running of PhenQuantConcAssigner - SKIP"
            )
            return

        logger.info(
            "'PhenotypeManager': completed quantitative phenotype data analysis for "
            "concentration data."
        )
