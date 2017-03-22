import pytest

@pytest.mark.xfail
def test_dupradar(sample1_se_dupradar):
    assert open(sample1_se_dupradar['dataframe']).readline().startswith('"ID"\t"geneLength"')
