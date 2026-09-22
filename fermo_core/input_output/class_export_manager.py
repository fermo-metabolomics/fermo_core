"""Organizes data export methods for dumping on disk.

Copyright (c) 2022-2024 Mitja Maximilian Zdouc, PhD

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

import json
import logging
import platform
from datetime import datetime
from typing import Any, Self

import networkx as nx
import pandas as pd
from pydantic import BaseModel

from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager
from fermo_core.input_output.class_summary_writer import SummaryWriter
from fermo_core.input_output.class_validation_manager import ValidationManager

logger = logging.getLogger("fermo_core")


class JsonExporter(BaseModel):
    """A Pydantic-based class for methods on exporting data in json format.

    Attributes:
        params: a Parameter object with specified user parameters
        stats: a Stats object containing general information
        features: a Repository object containing feature and general information
        samples: a Repository object containing sample information
        version: currently running version of fermo_core
        starttime: the date and time at start of fermo_core processing
        session: a dict for collecting data for json dump
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository
    version: str
    starttime: datetime
    session: dict = {}

    def export_metadata_json(self: Self):
        """Export metadata on analysis run."""
        self.session["metadata"] = {
            "fermo_core_version": self.version,
            "file_created_isoformat": datetime.now().isoformat(),
            "runtime_seconds": (datetime.now() - self.starttime).total_seconds(),
            "system": platform.system(),
            "version": platform.version(),
            "architecture": platform.architecture(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": platform.python_version(),
        }

    def export_params_json(self: Self):
        """Export data from params"""
        self.session["parameters"] = self.params.to_json()

    def export_stats_json(self: Self):
        """Export data from stats"""
        self.session["stats"] = self.stats.to_json()

    def export_features_json(self: Self):
        """Export data from features repository"""
        if self.stats.active_features is not None:
            self.session["general_features"] = {}
            for feature_id in self.stats.active_features:
                feature = self.features.get(feature_id)
                self.session["general_features"][feature_id] = feature.to_json()

    def export_samples_json(self: Self):
        """Export data from samples repository"""
        if self.stats.samples is not None:
            self.session["samples"] = {}
            for sample_id in self.stats.samples:
                sample = self.samples.get(sample_id)
                self.session["samples"][sample_id] = sample.to_json()

    def build_json_dict(self: Self):
        """Driver method to assemble data for json dump"""
        self.export_metadata_json()
        self.export_params_json()
        self.export_stats_json()
        self.export_features_json()
        self.export_samples_json()

    def return_session(self: Self) -> dict:
        """Return the generated session object to calling method"""
        return self.session


class CsvExporter(BaseModel):
    """A Pydantic-based class for methods on exporting data in csv format.

    Attributes:
        params: a Parameter object with specified user parameters
        stats: a Stats object containing general information
        features: a Repository object containing feature and general information
        samples: a Repository object containing sample information
        df: a Pandas dataframe acting as aggregator for data
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository
    df: Any

    def add_activity_info_csv(self: Self):
        """Iterate through active/inactive feature info and prepare for export"""

        def _add_activity_info(f_id: int) -> str:
            if f_id in self.stats.active_features:
                return "true"
            else:
                return "false"

        self.df["fermo:active"] = self.df["id"].map(
            lambda x: _add_activity_info(f_id=x)
        )

    def add_sample_info_csv(self: Self):
        """Iterate through feature sample information and prepare for export"""

        def _add_sample_info(f_id: int) -> str | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                return "|".join([sample for sample in feature.samples])
            return None

        def _add_sample_count_info(f_id: int) -> int | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                return len([sample for sample in feature.samples])
            return None

        self.df["fermo:samples"] = self.df["id"].map(lambda x: _add_sample_info(f_id=x))
        self.df["fermo:samples:count"] = self.df["id"].map(
            lambda x: _add_sample_count_info(f_id=x)
        )

    def add_blank_info_csv(self: Self):
        """Iterate through feature blank information and prepare for export"""

        def _add_blank(f_id: int) -> str | None:
            try:
                feature = self.features.get(f_id)
                if feature.blank is None:
                    return None
                elif feature.blank is True:
                    return "true"
                else:
                    return "false"
            except (TypeError, AttributeError, KeyError):
                return None

        self.df["fermo:isblank"] = self.df["id"].map(lambda x: _add_blank(f_id=x))

    def add_group_info_csv(self: Self):
        """Iterate through feature group information and prepare for export"""

        def _add_cat_groups(f_id: int, cat: str) -> str | None:
            try:
                feature = self.features.get(f_id)
                return "|".join(list(feature.groups[cat]))
            except (TypeError, AttributeError, KeyError):
                return None

        categories = [str(key) for key, val in self.stats.GroupMData.ctgrs.items()]
        if len(categories) != 0:
            for categ in categories:
                self.df[f"fermo:category:{categ}"] = self.df["id"].map(
                    lambda x: _add_cat_groups(f_id=x, cat=categ)
                )

    def add_networks_info_csv(self: Self):
        """Iterate through network information and prepare for export"""

        def _get_network_value(f_id: int, attr: str) -> str | None:
            try:
                attrs = attr.split(":")
                feature = self.features.get(f_id)
                return getattr(getattr(feature, attrs[1])[attrs[2]], attrs[3])
            except (TypeError, AttributeError, KeyError):
                return None

        if self.stats.networks is None or len(self.stats.networks) == 0:
            return

        networks = [f"fermo:networks:{nw}:network_id" for nw in self.stats.networks]
        for network in networks:
            self.df[network] = self.df["id"].map(
                lambda x: _get_network_value(x, network)
            )

    def add_adduct_info_csv(self: Self):
        """Iterate through adduct annotation information and add to df"""

        def _add_adduct_info(f_id: int) -> str | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                annot = []
                try:
                    for adct in feature.Annotations.adducts:
                        annot.append(
                            f"{adct.adduct_type}"
                            f"(partner_ID={adct.partner_id};"
                            f"adduct={adct.partner_adduct};"
                            f"m/z={adct.partner_mz})"
                        )
                    return "|".join(annot)
                except (TypeError, AttributeError, KeyError):
                    return None
            return None

        if not (
            self.params.AdductAnnotationParameters
            and self.params.AdductAnnotationParameters.activate_module
        ):
            return

        self.df["fermo:annotation:adducts"] = self.df["id"].map(
            lambda x: _add_adduct_info(f_id=x)
        )

    def add_loss_info_csv(self: Self):
        """Iterate through neutral loss annotation information and add to df"""

        def _add_loss_info(f_id: int) -> str | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                losses = []
                try:
                    for loss in feature.Annotations.losses:
                        losses.append(
                            f"'{loss.id}'"
                            f"(detected_loss={round(loss.loss_det, 4)};"
                            f"diff_ppm={round(loss.diff, 1)})"
                        )
                    return "|".join(losses)
                except (TypeError, AttributeError, KeyError):
                    return None
            return None

        if not (
            self.params.NeutralLossParameters
            and self.params.NeutralLossParameters.activate_module
        ):
            return

        self.df["fermo:annotation:neutral_losses"] = self.df["id"].map(
            lambda x: _add_loss_info(f_id=x)
        )

    def add_match_info_csv(self: Self):
        """Iterate through user library match annotation information and add to df"""

        def _add_match_info(f_id: int, var: str) -> str | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                matches = []
                try:
                    for match in feature.Annotations.matches:
                        if match.module == var:
                            matches.append(
                                f"'{match.id}'"
                                f"(score={match.score};"
                                f"algorithm={match.algorithm};"
                                f"diff_mz={match.diff_mz})"
                            )
                    return "|".join(matches)
                except (TypeError, AttributeError, KeyError):
                    return None
            return None

        modules = set()
        modules.add("MzmineAnnParser")
        if (
            self.params.SpectralLibMatchingDeepscoreParameters
            and self.params.SpectralLibMatchingDeepscoreParameters.activate_module
        ):
            modules.add("user_library_annotation")
        if (
            self.params.SpectralLibMatchingCosineParameters
            and self.params.SpectralLibMatchingCosineParameters.activate_module
        ):
            modules.add("user_library_annotation")
        if self.params.MS2QueryResultsParameters:
            modules.add("ms2query_annotation")
        if self.params.AsResultsParameters:
            modules.add("antismash_kcb_annotation")

        for module in modules:
            self.df[f"fermo:annotation:matches:{module}"] = self.df["id"].map(
                lambda x: _add_match_info(f_id=x, var=module)
            )

    def add_fragment_info_csv(self: Self):
        """Iterate through user library fragment annotation information and add to df"""

        def _add_fragment_info(f_id: int) -> str | None:
            if f_id in self.stats.active_features:
                feature = self.features.get(f_id)
                fragments = []
                try:
                    for frag in feature.Annotations.fragments:
                        fragments.append(
                            f"'{frag.id}'"
                            f"(detected_fragment={round(frag.frag_det, 4)};"
                            f"diff_ppm={round(frag.diff, 1)})"
                        )
                    return "|".join(fragments)
                except (TypeError, AttributeError, KeyError):
                    return None
            return None

        if not (
            self.params.FragmentAnnParameters
            and self.params.FragmentAnnParameters.activate_module
        ):
            return

        self.df["fermo:annotation:fragments"] = self.df["id"].map(
            lambda x: _add_fragment_info(f_id=x)
        )

    def add_phenotype_info_csv(self: Self):
        """Iterate through phenotype annotation and add to df"""

        if self.stats.phenotypes is None:
            return

        for categ in self.stats.phenotypes:
            self.df[f"fermo:phenotype:{categ.category}"] = self.df["id"].map(
                lambda x: "true" if x in categ.f_ids_positive else None
            )

    def build_csv_output(self: Self):
        """Assemble data for csv export"""
        self.add_activity_info_csv()
        self.add_sample_info_csv()
        self.add_blank_info_csv()
        self.add_group_info_csv()
        self.add_networks_info_csv()
        self.add_phenotype_info_csv()
        self.add_adduct_info_csv()
        self.add_loss_info_csv()
        self.add_match_info_csv()
        self.add_fragment_info_csv()

    def return_dfs(self: Self) -> tuple:
        """Return the generated df objects to calling method for export

        Returns:
            Tuple of a full dataframe and an abbreviated one
        """
        df_full = self.df.copy(deep=True)

        abbr_cols = ["id", "height", "area", "mz", "rt"]
        abbr_cols.extend([col for col in self.df if col.startswith("fermo")])
        df_abbr = self.df[abbr_cols].copy(deep=True)

        return df_full, df_abbr


class ExportManager(BaseModel):
    """A Pydantic-based class managing data export.

    Attributes:
        params: a Parameter object with specified user parameters
        stats: a Stats object containing general information
        features: a Repository object containing feature and general information
        samples: a Repository object containing sample information
    """

    params: ParameterManager
    stats: Stats
    features: Repository
    samples: Repository

    @staticmethod
    def log_start_module(file: str):
        """Log the start of the export of the corresponding file"""
        logger.debug(f"'ExportManager': started export of '{file}'.")

    @staticmethod
    def log_complete_module(file: str):
        """Log the completion of the export of the corresponding file"""
        logger.debug(f"'ExportManager': completed export of '{file}'.")

    def write_fermo_json(self: Self, version: str, starttime: datetime):
        """Write collected data in session dict into a json file on disk

        Arguments:
            version: a str indicating the currently running version of fermo_core
            starttime: the date and time at start of fermo_core processing
        """
        self.log_start_module("fermo.session.json")

        json_exporter = JsonExporter(
            params=self.params,
            stats=self.stats,
            features=self.features,
            samples=self.samples,
            version=version,
            starttime=starttime,
        )
        json_exporter.build_json_dict()
        session = json_exporter.return_session()

        session_path = self.params.OutputParameters.directory_path.joinpath(
            "out.fermo.session.json"
        )

        with open(session_path, "w", encoding="utf-8") as outfile:
            outfile.write(json.dumps(session, indent=2, ensure_ascii=False))

        ValidationManager().validate_output_created(session_path)

        self.log_complete_module("fermo.session.json")

    def write_csv_output(self: Self):
        """Write modified peaktable as csv on disk"""
        self.log_start_module("fermo.full.csv/fermo.abbrev.csv")

        csv_exporter = CsvExporter(
            params=self.params,
            stats=self.stats,
            features=self.features,
            samples=self.samples,
            df=pd.read_csv(self.params.PeaktableParameters.filepath),
        )
        csv_exporter.build_csv_output()
        df_full, df_abbr = csv_exporter.return_dfs()

        path_df_full = self.params.OutputParameters.directory_path.joinpath(
            "out.fermo.full.csv"
        )
        path_df_abbr = self.params.OutputParameters.directory_path.joinpath(
            "out.fermo.abbrev.csv"
        )

        df_full.to_csv(path_df_full, encoding="utf-8", index=False, sep=",")
        df_abbr.to_csv(path_df_abbr, encoding="utf-8", index=False, sep=",")

        ValidationManager().validate_output_created(path_df_full)
        ValidationManager().validate_output_created(path_df_abbr)

        self.log_complete_module("fermo.full.csv/fermo.abbrev.csv")

    def write_cytoscape_output(self: Self):
        """Write Cytoscape-compatible graphml output if networking was performed"""

        self.log_start_module(".graphml")

        if self.stats.networks is None:
            logger.warning(
                "'ExportManager': no spectral similarity networks generated - SKIP"
            )
            return

        for network in self.stats.networks:
            path_graphml = self.params.OutputParameters.directory_path.joinpath(
                f"out.fermo.{network}.graphml"
            )
            nx.write_graphml(self.stats.networks[network].network, path_graphml)
            ValidationManager().validate_output_created(path_graphml)

        self.log_complete_module(".graphml")

    def write_summary_output(self: Self):
        """Write a human-readable summary of steps performed."""
        dst = self.params.OutputParameters.directory_path.joinpath(
            "out.fermo.summary.txt"
        )
        self.log_start_module("summary.txt")
        summary_writer = SummaryWriter(params=self.params, destination=dst)
        summary_writer.assemble_summary()
        summary_writer.write_summary()
        ValidationManager().validate_output_created(dst)
        self.log_complete_module("summary.txt")

    def run(self: Self, version: str, starttime: datetime):
        """Call export methods based on user-input

        Arguments:
            version: a str indicating the currently running version of fermo_core
            starttime: the date and time at start of fermo_core processing
        """
        self.write_fermo_json(version, starttime)
        self.write_csv_output()
        self.write_cytoscape_output()
        self.write_summary_output()
