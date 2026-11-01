"""ISOTXS assembly and storage, without running OpenMC."""

import numpy as np
from armi.nuclearDataIO.cccc import isotxs

from armi_openmc import mgxs
from armi_openmc.materials import NuclideMapper
from armi_openmc.tests.conftest import crossSections, needsData


def test_bandStructureWithUpscatter():
    # [to, from]: group 0 fed by itself; group 1 by 0 and 1 and (upscatter) by 2; group 2 by 1, 2.
    s = np.array([[1.0, 0.0, 0.0], [0.5, 1.0, 0.2], [0.0, 0.3, 1.0]])
    jband, jj, up, down = mgxs._bandStructure(s)
    assert jband == [1, 3, 2]
    assert jj == [1, 2, 1]
    assert (up, down) == (1, 1)


class _Model:
    xsID = "AA"
    temperatureK = 600.0


def _constants(ng):
    rng = np.random.default_rng(1)
    scatter = np.tril(rng.uniform(0.1, 1.0, (ng, ng)))
    scatter[0, 1] = 0.05  # one upscatter term
    return mgxs.GroupConstants(
        total=rng.uniform(5, 10, ng),
        transport=rng.uniform(4, 9, ng),
        fission=rng.uniform(1, 2, ng),
        nu=np.full(ng, 2.4),
        chi=np.linspace(1, 0, ng) / np.linspace(1, 0, ng).sum(),
        capture=rng.uniform(0.1, 1, ng),
        extras={k: np.zeros(ng) for k in mgxs.ISOTXS_EXTRAS.values()},
        scatter=scatter,
        efiss=3.2e-11,
        density=0.01,
    )


@needsData
def test_libraryRoundTrip(tmp_path):
    ng = 33
    gc = _constants(ng)
    info = {"velocity": np.linspace(1e9, 1e5, ng), "flux": np.ones(ng), "keff": None}
    lib = mgxs.buildLibrary(
        [(_Model(), {"U235": gc}, info)], "ANL33", NuclideMapper(crossSections())
    )
    path = tmp_path / "ISOTXS"
    mgxs.writeLibrary(lib, str(path))
    back = isotxs.readBinary(str(path))
    nuc = back["U235AA"]
    assert np.allclose(nuc.micros.fission, gc.fission)
    assert np.allclose(nuc.micros.chi, gc.chi)
    assert np.allclose(nuc.micros.elasticScatter.toarray(), gc.scatter)
    assert back.isotxsMetadata["maxUpScatterGroups"] == 1
    assert nuc.isotxsMetadata["ecapt"] > 0  # from the capture Q value in the HDF5 data
