import pytest

from fermo_core.data_analysis.phenotype_manager.class_phen_qual_assigner import (
    PhenQualAssigner,
)
from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Feature,
    SampleInfo,
)
from fermo_core.data_processing.builder_sample.dataclass_sample import Sample
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import PhenoData, SamplePhenotype, Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager
from fermo_core.input_output.param_handlers import (
    PhenoQualAssgnParams,
    PhenotypeParameters,
)


@pytest.fixture
def phen_qual():
    phen_qual = PhenQualAssigner(
        params=ParameterManager(),
        stats=Stats(),
        features=Repository(),
        samples=Repository(),
    )
    f1 = Feature(
        f_id=1,
        samples={
            "s1",
        },
        area_per_sample=[
            SampleInfo(s_id="s1", value=10),
        ],
        height_per_sample=[
            SampleInfo(s_id="s1", value=100),
        ],
    )
    f2 = Feature(
        f_id=2,
        samples={"s1", "s2"},
        area_per_sample=[
            SampleInfo(s_id="s1", value=10),
            SampleInfo(s_id="s2", value=1),
        ],
        height_per_sample=[
            SampleInfo(s_id="s1", value=100),
            SampleInfo(s_id="s2", value=10),
        ],
    )
    f3 = Feature(
        f_id=3,
        samples={
            "s2",
        },
        area_per_sample=[
            SampleInfo(s_id="s2", value=10),
        ],
        height_per_sample=[
            SampleInfo(s_id="s2", value=100),
        ],
    )
    phen_qual.features.add(1, f1)
    phen_qual.features.add(2, f2)
    phen_qual.features.add(3, f3)
    s1 = Sample(s_id="s1", feature_ids={1, 2})
    s2 = Sample(s_id="s2", feature_ids={2, 3})
    phen_qual.samples.add("s1", s1)
    phen_qual.samples.add("s2", s2)
    phen_qual.stats.active_features = {1, 2, 3}
    phen_qual.stats.samples = ("s1", "s2")
    phen_qual.stats.phenotypes = [
        PhenoData(
            datatype="qualitative",
            category="qualitative",
            s_negative={
                "s2",
            },
            s_phen_data=[SamplePhenotype(s_id="s1")],
        )
    ]
    phen_qual.params.PhenoQualAssgnParams = PhenoQualAssgnParams(
        **{
            "activate_module": True,
            "factor": 5,
            "algorithm": "minmax",
            "value": "area",
            "p_val_cutoff": 0.0,
        }
    )
    phen_qual.params.PhenotypeParameters = PhenotypeParameters(
        **{
            "filepath": "tests/test_data/test.bioactivity.qualitative.csv",
            "format": "qualitative",
        }
    )
    return phen_qual


def test_collect_sets_valid(phen_qual):
    phen_qual.collect_sets()
    assert len(phen_qual.f_ids_intersect) == 1


def test_collect_sets_invalid(phen_qual):
    phen_qual.samples.entries["s1"].feature_ids = {
        1,
    }
    phen_qual.collect_sets()
    assert len(phen_qual.f_ids_intersect) == 0


def test_get_value_valid(phen_qual):
    vals = phen_qual.get_value(2, {"s1", "s2"})
    assert len(vals) == 2


def test_get_value_invalid(phen_qual):
    with pytest.raises(RuntimeError):
        phen_qual.get_value(1, {"s2"})


def test_bin_intersection_positive(phen_qual):
    phen_qual.f_ids_intersect = {
        2,
    }
    phen_qual.bin_intersection()
    assert len(phen_qual.stats.phenotypes[0].f_ids_positive) == 1


def test_bin_intersection_negative(phen_qual):
    phen_qual.f_ids_intersect = {
        2,
    }
    phen_qual.features.entries[2].area_per_sample = [
        SampleInfo(s_id="s1", value=10),
        SampleInfo(s_id="s2", value=10),
    ]
    phen_qual.bin_intersection()
    assert len(phen_qual.stats.phenotypes[0].f_ids_positive) == 0


def test_determine_active(phen_qual):
    phen_qual.params.PhenoQualAssgnParams.p_val_cutoff = 0.05
    feature = phen_qual.determine_active(
        feature=phen_qual.features.get(1),
        fct=10,
        act=[66, 77, 88, 99, 111],
        inact=[6, 7, 8, 9],
    )
    assert feature.Annotations.phenotypes[0].p_value is not None
