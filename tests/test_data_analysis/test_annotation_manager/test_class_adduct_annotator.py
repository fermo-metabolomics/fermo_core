import networkx as nx
import pytest

from fermo_core.data_analysis.annotation_manager.class_adduct_annotator import (
    AdductAnnotator,
)
from fermo_core.data_processing.builder_feature.dataclass_feature import (
    Adduct,
    Feature,
    SimNetworks,
)
from fermo_core.data_processing.class_repository import Repository
from fermo_core.data_processing.class_stats import Stats
from fermo_core.input_output.class_parameter_manager import ParameterManager
from fermo_core.input_output.param_handlers import AdductAnnotationParameters


@pytest.fixture
def adduct_annotator(
    parameter_instance, stats_instance, feature_instance, sample_instance
):
    return AdductAnnotator(
        params=parameter_instance,
        stats=stats_instance,
        features=feature_instance,
        samples=sample_instance,
    )


@pytest.fixture
def adduct_annotator_min():
    params = ParameterManager(
        AdductAnnotationParameters=AdductAnnotationParameters(
            **{"activate_module": True, "mass_dev_ppm": 10.0}
        )
    )
    features = Repository()
    feature1 = Feature(f_id=1, mz=100.0)
    features.add(1, feature1)
    feature2 = Feature(f_id=2, mz=100.0)
    features.add(2, feature2)
    stats = Stats(active_features={1, 2})
    return AdductAnnotator(
        params=params, stats=stats, features=features, samples=Repository()
    )


def test_return_features_valid(adduct_annotator_min):
    assert adduct_annotator_min.return_features() is not None


def test_run_analysis_valid(adduct_annotator):
    adduct_annotator.run_analysis()
    features, stats = adduct_annotator.return_features()
    assert features.entries[131].Annotations.adducts[0] is not None


def test_create_network_singletons_valid(adduct_annotator_min):
    """Condiction: feature have no adducts - singletons in graph"""
    adduct_annotator_min.create_network()
    assert (
        len(
            list(
                nx.connected_components(
                    adduct_annotator_min.stats.networks["ion_identity"].network
                )
            )
        )
        == 2
    )


def test_create_network_subgraph_valid(adduct_annotator_min):
    """Condition: features are related and in same subgraph"""
    adduct_annotator_min.features.entries[1].mz = 415.2098
    adduct_annotator_min.features.entries[2].mz = 437.1912
    adduct_annotator_min.sodium_adduct(1, 2, "sample1")
    adduct_annotator_min.create_network()
    assert (
        len(
            list(
                nx.connected_components(
                    adduct_annotator_min.stats.networks["ion_identity"].network
                )
            )
        )
        == 1
    )


def test_create_network_two_subgraphs(adduct_annotator_min):
    """Condition: two subgraphs of two features each"""
    adduct_annotator_min.features.entries[1].mz = 415.2098
    adduct_annotator_min.features.entries[2].mz = 437.1912
    adduct_annotator_min.sodium_adduct(1, 2, "sample1")

    adduct_annotator_min.features.add(3, Feature(f_id=3, mz=1648.4547))
    adduct_annotator_min.features.add(4, Feature(f_id=4, mz=1649.4578))
    adduct_annotator_min.stats.active_features.update([3, 4])
    adduct_annotator_min.plus1_isotope(3, 4, "sample1")
    adduct_annotator_min.create_network()
    assert (
        len(
            list(
                nx.connected_components(
                    adduct_annotator_min.stats.networks["ion_identity"].network
                )
            )
        )
        == 2
    )
    assert adduct_annotator_min.stats.networks["ion_identity"].summary == {
        0: {1, 2},
        1: {3, 4},
    }
    assert isinstance(
        adduct_annotator_min.features.entries[1].networks["ion_identity"], SimNetworks
    )


def test_add_adduct_info_valid(adduct_annotator_min):
    feature = adduct_annotator_min.add_adduct_info(Feature())
    assert feature.Annotations is not None


def test_annotate_spec_features_valid(adduct_annotator):
    adduct_annotator.annotate_adducts_pos("5458_5457_mod.mzXML")
    features, stats = adduct_annotator.return_features()
    assert features.entries[131].Annotations.adducts[0] is not None


def test_dereplicate_adducts_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 415.2098
    adduct_annotator_min.features.entries[2].mz = 437.1912
    adduct_annotator_min.sodium_adduct(1, 2, "sample1")
    adduct_annotator_min.features.entries[1].Annotations.adducts.append(
        Adduct(
            adduct_type="[M+H]+",
            partner_adduct="[M+Na]+",
            partner_id=2,
            partner_mz=437.1912,
            diff_ppm=12.0,
            sample="sample2",
        )
    )
    adduct_annotator_min.dereplicate_adducts()
    features, stats = adduct_annotator_min.return_features()
    adduct_dict = features.entries[1].Annotations.adducts[0].to_json()
    assert len(adduct_dict["samples"]) == 2


def test_sodium_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 415.2098
    adduct_annotator_min.features.entries[2].mz = 437.1912
    assert adduct_annotator_min.sodium_adduct(1, 2, "sample1")


def test_dimer_sodium_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 415.2098
    adduct_annotator_min.features.entries[2].mz = 851.39487
    assert adduct_annotator_min.dimer_sodium_adduct(1, 2, "sample1")


def test_triple_h_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1510.4198
    adduct_annotator_min.features.entries[2].mz = 504.1447
    assert adduct_annotator_min.triple_h_adduct(1, 2, "sample1")


def test_plus1_isotope_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 1649.4578
    assert adduct_annotator_min.plus1_isotope(1, 2, "sample1")


def test_plus2_isotope_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 1650.4653
    assert adduct_annotator_min.plus2_isotope(1, 2, "sample1")


def test_plus3_isotope_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 1651.4547
    assert adduct_annotator_min.plus3_isotope(1, 2, "sample1")


def test_plus4_isotope_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 1652.4539
    assert adduct_annotator_min.plus4_isotope(1, 2, "sample1")


def test_plus5_isotope_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 1653.4754
    assert adduct_annotator_min.plus5_isotope(1, 2, "sample1")


def test_double_plus1_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 825.2326
    assert adduct_annotator_min.double_plus1(1, 2, "sample1")


def test_double_plus2_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 825.7343
    assert adduct_annotator_min.double_plus2(1, 2, "sample1")


def test_double_plus3_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 826.2360
    assert adduct_annotator_min.double_plus3(1, 2, "sample1")


def test_double_plus4_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 826.7377
    assert adduct_annotator_min.double_plus4(1, 2, "sample1")


def test_double_plus5_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.4547
    adduct_annotator_min.features.entries[2].mz = 827.2393
    assert adduct_annotator_min.double_plus5(1, 2, "sample1")


def test_iron56_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 843.4772
    adduct_annotator_min.features.entries[2].mz = 896.3883
    assert adduct_annotator_min.iron56(1, 2, "sample1")


def test_dimer_double_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1510.4198
    adduct_annotator_min.features.entries[2].mz = 755.7153
    assert adduct_annotator_min.dimer_double(1, 2, "sample1")


def test_ammonium_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 409.29477
    adduct_annotator_min.features.entries[2].mz = 426.321
    assert adduct_annotator_min.ammonium(1, 2, "sample1")


def test_potassium_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 409.29477
    adduct_annotator_min.features.entries[2].mz = 447.251
    assert adduct_annotator_min.potassium(1, 2, "sample1")


def test_water_add_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 409.29477
    adduct_annotator_min.features.entries[2].mz = 427.30588
    assert adduct_annotator_min.water_add(1, 2, "sample1")


def test_water_loss_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 409.29477
    adduct_annotator_min.features.entries[2].mz = 391.284
    assert adduct_annotator_min.water_loss(1, 2, "sample1")


def test_chloride_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 852.323614
    adduct_annotator_min.features.entries[2].mz = 888.300292
    assert adduct_annotator_min.chloride_adduct(1, 2, "sample1")


def test_double_dimer_pair_neg_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1648.47
    adduct_annotator_min.features.entries[2].mz = 823.73
    assert adduct_annotator_min.double_dimer_pair_neg(1, 2, "sample1")


def test_bicarbonate_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 852.323614
    adduct_annotator_min.features.entries[2].mz = 914.32345
    assert adduct_annotator_min.bicarbonate_adduct(1, 2, "sample1")


def test_tfa_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 852.323614
    adduct_annotator_min.features.entries[2].mz = 966.316476
    assert adduct_annotator_min.tfa_adduct(1, 2, "sample1")


def test_acetate_adduct_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 852.323614
    adduct_annotator_min.features.entries[2].mz = 912.344741
    assert adduct_annotator_min.acetate_adduct(1, 2, "sample1")


def test_double_quadruple_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1405.58893
    adduct_annotator_min.features.entries[2].mz = 703.29907
    assert adduct_annotator_min.double_quadruple(1, 2, "sample1")


def test_double_triple_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 1405.58893
    adduct_annotator_min.features.entries[2].mz = 937.39628
    assert adduct_annotator_min.double_triple(1, 2, "sample1")


def test_quadruple_triple_valid(adduct_annotator_min):
    adduct_annotator_min.features.entries[1].mz = 703.29907
    adduct_annotator_min.features.entries[2].mz = 937.39628
    assert adduct_annotator_min.quadruple_triple(1, 2, "sample1")
