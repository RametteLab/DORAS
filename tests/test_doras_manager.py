import pytest
from src.doras_manager import DorasManager
from src.config import DorasParams, Phase


@pytest.fixture
def mock_params():
    return DorasParams(
        genome_size=5000000,
        selected_db="test_db",
        sample_names=["BC01", "BC02"],
        mlst_genes_path="test_mlst_path",
        fastq_files_path="test_fastq_path",
    )


def test_doras_manager_initialization(mock_params, tmp_path, monkeypatch):
    # Using monkeypatch instead of patch
    class MockQuery:
        def __init__(self, **kwargs):
            pass

        def initialize_extended_ref(self):
            pass

    monkeypatch.setattr("src.doras_manager.ExtensionQuery", MockQuery)

    manager = DorasManager(
        barcodes=["BC01", "BC02"],
        folder_path=str(tmp_path),
        mlst_ref_genes_path="test_mlst_path",
        mapq=60,
        test_mode=True,
        overwrite=False,
        clean_up=False,
        target_extension_size=3000,
        params=mock_params,
    )

    assert manager.barcodes == ["BC01", "BC02"]
    assert manager.output_dir == tmp_path / "DORAS"
    assert manager.params == mock_params
    assert len(manager.processors) == 2


def test_setup_processors(mock_params, tmp_path, monkeypatch):
    # Mocking ExtensionQuery to track if it's called
    class MockQuery:
        called = False

        def __init__(self, **kwargs):
            MockQuery.called = True

        def initialize_extended_ref(self):
            pass

    monkeypatch.setattr("src.doras_manager.ExtensionQuery", MockQuery)

    manager = DorasManager(
        barcodes=["BC01"],
        folder_path=str(tmp_path),
        mlst_ref_genes_path="test_mlst_path",
        mapq=60,
        test_mode=False,
        overwrite=False,
        clean_up=False,
        target_extension_size=3000,
        params=mock_params,
    )

    # Processors are already setup in __init__
    assert len(manager.processors) == 1
    assert MockQuery.called


def test_doras_manager_phase_initialization(mock_params, tmp_path, monkeypatch):
    class MockQuery:
        def __init__(self, **kwargs):
            pass

        def initialize_extended_ref(self):
            pass

    monkeypatch.setattr("src.doras_manager.ExtensionQuery", MockQuery)

    manager = DorasManager(
        barcodes=["BC01"],
        folder_path=str(tmp_path),
        mlst_ref_genes_path="test_mlst_path",
        mapq=60,
        test_mode=True,
        overwrite=False,
        clean_up=False,
        target_extension_size=3000,
        params=mock_params,
        phase=Phase.QUERY_ST,
    )
    assert manager.phase == Phase.QUERY_ST
