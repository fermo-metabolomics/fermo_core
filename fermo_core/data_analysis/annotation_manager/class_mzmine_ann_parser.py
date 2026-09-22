"""Retrieves optional annotations from mzmine peaktable

Copyright (c) 2022-present Mitja Maximilian Zdouc, PhD

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

import pandas as pd
from pydantic import BaseModel

from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Adduct,
    Annotations,
    Feature,
    Match,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager

logger = logging.getLogger("fermo_core")


class MzmineAnnParser(BaseModel):
    """Retrieves optional annotations from mzmine peaktable file

    Attributes:
        params: ParameterManager object, holds user-provided parameters
        stats: Stats object, holds info on currently active features (not filtered out)
        features: Repository object, holds "General Feature" objects
        accepted: parsed columns
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    f: Feature = Feature()
    accepted: tuple = (
        "ion_identities:ion_identities",
        "spectral_db_matches:spectral_db_matches",
    )

    def return_attributes(self) -> Repository:
        """Returns modified attributes

        Returns:
            Modified Feature object
        """
        return self.features

    def run(self) -> None:
        """Orchestrate functions"""
        logger.info("'MzmineAnnParser': Start analysis")

        if not self.params.PeaktableParameters.format in ["mzmine3", "mzmine4"]:
            logger.info("'MzmineAnnParser': not a mzmine table - SKIP")
            return

        df = pd.read_csv(self.params.PeaktableParameters.filepath)

        if not any(col in df.columns for col in self.accepted):
            logger.info(
                f"'MzmineAnnParser': did not detect optional annotation columns '{self.accepted}' - SKIP"
            )
            return

        for _, row in df.iterrows():
            if not self.contains_values(row):
                continue

            if not int(row["id"]) in self.stats.active_features:
                logger.warning(
                    f"'MzmineAnnParser': {int(row['id'])} not found in 'active_features' - SKIP"
                )
                continue

            self.f = self.features.get(int(row["id"]))

            if not self.f.Annotations:
                self.f.Annotations = Annotations()

            self.parse_ion_ident(row, df)
            self.parse_spectral_db(row)

            self.features.modify(int(row["id"]), self.f)

        logger.info("'MzmineAnnParser': Complete analysis")

    def contains_values(self, r: pd.Series) -> bool:
        """Check if any value is not NaN

        Args:
            r: a mzmine peaktable row
        """
        return any(pd.notnull(r.get(i)) for i in self.accepted)

    def parse_ion_ident(self, r: pd.Series, df: pd.DataFrame) -> None:
        """Retrieve ion identity information

        Args:
            r: a mzmine peaktable row
            df: the full dataframe for reference
        """
        if not pd.notnull(r.get("ion_identities:ion_identities")):
            return

        if not self.f.Annotations.adducts:
            self.f.Annotations.adducts = []

        try:
            for partner_id in r["ion_identities:partner_row_ids"].split(";"):
                p_id = int(partner_id)
                if int(r["id"]) == p_id:
                    continue

                self.f.Annotations.adducts.append(
                    Adduct(
                        adduct_type=f'{r["ion_identities:ion_identities"]}(mzmine)',
                        partner_id=p_id,
                        partner_adduct=f"{df.loc[df['id'] == p_id, 'ion_identities:ion_identities'].values[0]}(mzmine)",
                        partner_mz=df.loc[df["id"] == p_id, "mz"].values[0],
                        diff_ppm="0",
                    )
                )
        except Exception as e:
            logger.warning(
                f"MzmineAnnParser: could not find ion identity annotation: {e!s}"
            )

    def parse_spectral_db(self, r: pd.Series) -> None:
        """Retrieve spectral db annotation information

        Args:
            r: a mzmine peaktable row
        """

        def _isnan(val) -> str:
            if pd.isna(val):
                return "N/A"
            else:
                return val

        if not pd.notnull(r.get("spectral_db_matches:spectral_db_matches")):
            return

        if not self.f.Annotations.matches:
            self.f.Annotations.matches = []

        try:
            self.f.Annotations.matches.append(
                Match(
                    id=_isnan(r.get("spectral_db_matches:compound_name", "N/A")),
                    library=self.params.PeaktableParameters.filepath.name,
                    algorithm="mzmine",
                    score=_isnan(r.get("spectral_db_matches:similarity_score", "N/A")),
                    mz=_isnan(r.get("spectral_db_matches:precursor_mz", "N/A")),
                    diff_mz=_isnan(r.get("spectral_db_matches:mz_diff", "N/A")),
                    module="MzmineAnnParser",
                    smiles=_isnan(r.get("spectral_db_matches:smiles", "N/A")),
                    inchikey=_isnan(r.get("spectral_db_matches:inchi", "N/A")),
                )
            )
        except IndexError:
            logger.warning(
                f"MzmineAnnParser: could not assign db annotation for feature {r['id']}"
            )
