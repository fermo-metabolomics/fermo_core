import numpy as np
import pytest

from fermo_core.data_analysis.phenotype_manager.class_phen_quant_ass import PhenQuantAss
from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Feature,
    SampleInfo,
)
from fermo_core.data_processing.builder_sample.dataclass_sample import Sample
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import PhenoData, SamplePhenotype, Stats


@pytest.fixture
def phen_quant_perc():
    phen_quant_perc = PhenQuantAss(
        coeff_cutoff=0.7,
        p_val_cutoff=0.05,
        fdr_corr="bonferroni",
        alg="pearson",
        mode="percentage",
        stats=Stats(),
        features=Repository(),
        samples=Repository(),
    )
    f1 = Feature(
        f_id=1,
        samples={"s1", "s2", "s3", "s4"},
        area_per_sample=[
            SampleInfo(s_id="s1", value=100),
            SampleInfo(s_id="s2", value=10),
            SampleInfo(s_id="s3", value=1000),
            SampleInfo(s_id="s4", value=50),
        ],
    )
    f2 = Feature(
        f_id=2,
        samples={"s1", "s2", "s3", "s4"},
        area_per_sample=[
            SampleInfo(s_id="s1", value=10),
            SampleInfo(s_id="s2", value=100),
            SampleInfo(s_id="s3", value=1000),
            SampleInfo(s_id="s4", value=10000),
        ],
    )
    phen_quant_perc.features.add(1, f1)
    phen_quant_perc.features.add(2, f2)
    s1 = Sample(s_id="s1", feature_ids={1, 2})
    s2 = Sample(s_id="s2", feature_ids={1, 2})
    s3 = Sample(s_id="s3", feature_ids={1, 2})
    s4 = Sample(s_id="s4", feature_ids={1, 2})
    phen_quant_perc.samples.add("s1", s1)
    phen_quant_perc.samples.add("s2", s2)
    phen_quant_perc.samples.add("s3", s3)
    phen_quant_perc.samples.add("s4", s4)
    phen_quant_perc.stats.active_features = {1, 2}
    phen_quant_perc.stats.samples = ("s1", "s2", "s3", "s4")
    phen_quant_perc.stats.phenotypes = [
        PhenoData(
            datatype="quantitative-percentage",
            category="assay:assay1",
            s_phen_data=[
                SamplePhenotype(s_id="s1", value=0),
                SamplePhenotype(s_id="s2", value=1),
                SamplePhenotype(s_id="s3", value=10),
                SamplePhenotype(s_id="s4", value=100),
            ],
        ),
        PhenoData(
            datatype="quantitative-percentage",
            category="assay:assay2",
            s_phen_data=[
                SamplePhenotype(s_id="s1", value=100),
                SamplePhenotype(s_id="s2", value=10),
                SamplePhenotype(s_id="s3", value=1000),
                SamplePhenotype(s_id="s4", value=50),
            ],
        ),
    ]
    return phen_quant_perc


@pytest.fixture
def phen_quant_conc():
    phen_quant_conc = PhenQuantAss(
        coeff_cutoff=0.7,
        p_val_cutoff=0.05,
        fdr_corr="bonferroni",
        alg="pearson",
        mode="concentration",
        stats=Stats(),
        features=Repository(),
        samples=Repository(),
    )
    f1 = Feature(
        f_id=1,
        samples={"s1", "s2", "s3", "s4", "s5"},
        area_per_sample=[
            SampleInfo(s_id="s1", value=100),
            SampleInfo(s_id="s2", value=10),
            SampleInfo(s_id="s3", value=1000),
            SampleInfo(s_id="s4", value=50),
            SampleInfo(s_id="s5", value=300),
        ],
    )
    f2 = Feature(
        f_id=2,
        samples={"s1", "s2", "s3", "s4", "s5"},
        area_per_sample=[
            SampleInfo(s_id="s1", value=0),
            SampleInfo(s_id="s2", value=20),
            SampleInfo(s_id="s3", value=40),
            SampleInfo(s_id="s4", value=80),
            SampleInfo(s_id="s5", value=160),
        ],
    )
    phen_quant_conc.features.add(1, f1)
    phen_quant_conc.features.add(2, f2)
    phen_quant_conc.samples.add("s1", Sample(s_id="s1", feature_ids={1, 2}))
    phen_quant_conc.samples.add("s2", Sample(s_id="s2", feature_ids={1, 2}))
    phen_quant_conc.samples.add("s3", Sample(s_id="s3", feature_ids={1, 2}))
    phen_quant_conc.samples.add("s4", Sample(s_id="s4", feature_ids={1, 2}))
    phen_quant_conc.samples.add("s5", Sample(s_id="s5", feature_ids={1, 2}))
    phen_quant_conc.stats.active_features = (1, 2)
    phen_quant_conc.stats.samples = ("s1", "s2", "s3", "s4", "s5")
    phen_quant_conc.stats.phenotypes = [
        PhenoData(
            datatype="quantitative-concentration",
            category="assay:assay1",
            s_phen_data=[
                SamplePhenotype(s_id="s1", value=0),
                SamplePhenotype(s_id="s2", value=256),
                SamplePhenotype(s_id="s3", value=128),
                SamplePhenotype(s_id="s4", value=32),
                SamplePhenotype(s_id="s5", value=8),
            ],
        )
    ]
    return phen_quant_conc


def test_find_relevant_f_ids_perc(phen_quant_perc):
    phen_quant_perc.find_relevant_f_ids()
    assert len(phen_quant_perc.relevant_f_ids) == 2


def test_calculate_correlation_valid_perc(phen_quant_perc):
    phen_quant_perc.find_relevant_f_ids()
    phen_quant_perc.calculate_correlation()
    assert round(phen_quant_perc.assays["assay:assay1"][1]["corr"], 3) == -0.249


def test_correct_p_val_perc(phen_quant_perc):
    phen_quant_perc.find_relevant_f_ids()
    phen_quant_perc.calculate_correlation()
    phen_quant_perc.correct_p_val()
    assert phen_quant_perc.assays["assay:assay1"][1]["adjusted_p"] == 1.0


def test_assign_results_perc(phen_quant_perc):
    phen_quant_perc.find_relevant_f_ids()
    phen_quant_perc.calculate_correlation()
    phen_quant_perc.correct_p_val()
    phen_quant_perc.assign_results()
    assert (
        phen_quant_perc.features.entries[2].Annotations.phenotypes[0].category
        == "assay:assay1"
    )


def test_run_analysis_valid_perc(phen_quant_perc):
    phen_quant_perc.run_analysis()
    assert phen_quant_perc.features.entries[2].Annotations.phenotypes[0] is not None


def test_run_analysis_invalid_perc(phen_quant_perc):
    phen_quant_perc.stats.phenotypes = None
    with pytest.raises(RuntimeError):
        phen_quant_perc.run_analysis()


def test_find_relevant_f_ids_conc(phen_quant_conc):
    phen_quant_conc.find_relevant_f_ids()
    assert len(phen_quant_conc.relevant_f_ids) == 2


def test_calculate_correlation_valid_conc(phen_quant_conc):
    phen_quant_conc.find_relevant_f_ids()
    phen_quant_conc.calculate_correlation()
    assert round(phen_quant_conc.assays["assay:assay1"][1]["corr"], 3) == -0.025


def test_correct_p_val_conc(phen_quant_conc):
    phen_quant_conc.find_relevant_f_ids()
    phen_quant_conc.calculate_correlation()
    phen_quant_conc.correct_p_val()
    assert phen_quant_conc.assays["assay:assay1"][1]["adjusted_p"] == 1.0


def test_assign_results_conc(phen_quant_conc):
    phen_quant_conc.find_relevant_f_ids()
    phen_quant_conc.calculate_correlation()
    phen_quant_conc.correct_p_val()
    phen_quant_conc.assign_results()
    assert (
        phen_quant_conc.features.entries[2].Annotations.phenotypes[0].category
        == "assay:assay1"
    )


def test_run_analysis_valid_conc(phen_quant_conc):
    phen_quant_conc.run_analysis()
    assert phen_quant_conc.features.entries[2].Annotations.phenotypes[0] is not None


def test_run_analysis_invalid_conc(phen_quant_conc):
    phen_quant_conc.stats.phenotypes = None
    with pytest.raises(RuntimeError):
        phen_quant_conc.run_analysis()


def test_valid_constant_valid(phen_quant_conc):
    assert phen_quant_conc.valid_constant([0, 1, 2, 3, 4])


def test_valid_constant_invalid(phen_quant_conc):
    assert phen_quant_conc.valid_constant([1, 1, 1, 1, 1]) is False
    assert phen_quant_conc.valid_constant([1, np.nan, 1, 1, 1]) is False


def test_pearson_percentage_valid(phen_quant_perc):
    x, y = phen_quant_perc.pearson([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], mode="percentage")
    assert x == 1


def test_pearson_concentration_valid(phen_quant_conc):
    x, y = phen_quant_conc.pearson(
        [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], mode="concentration"
    )
    assert round(x, 2) == 0.9


def test_spearman_percentage_valid(phen_quant_perc):
    x, y = phen_quant_perc.pearson([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], mode="percentage")
    assert x == 1


def test_spearman_concentration_valid(phen_quant_conc):
    x, y = phen_quant_conc.pearson(
        [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], mode="concentration"
    )
    assert round(x, 2) == 0.9


def test_spearman_perm_percentage_valid(phen_quant_perc):
    x, y = phen_quant_perc.spearman_permutation(
        [1, 2, 3, 4, 5], [1, 2, 3, 4, 5], mode="percentage"
    )
    assert round(x, 2) == 1


def test_spearman_perm_concentration_valid(phen_quant_conc):
    x, y = phen_quant_conc.spearman_permutation(
        [1, 2, 3, 4, 5], [5, 4, 3, 2, 1], mode="concentration"
    )
    assert round(x, 2) == 1.0
